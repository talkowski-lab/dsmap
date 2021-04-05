#######################
#    DSMap Project    #
#######################
#
# TrainMuModel.wdl
#
# Intersect SVs with bin-pairs and train mutation rate model (with optional evaluation)
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"
import "ShardAndPredictMu.wdl"


workflow TrainMuModel {
  input {
    # General inputs
    String vcf
    File vcf_idx
    String pairs_bucket
    String? pairs_bed_prefix
    File contigs_fai
    File training_mask
    String model
    Int shard_size_apply_mu
    String cnv
    String prefix

    # Diagnostics options
    Boolean run_diagnostics = true

    # Dockers
    String athena_docker
    String athena_cloud_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_intersect_svs
    RuntimeAttr? runtime_attr_apply_training_mask
    RuntimeAttr? runtime_attr_diagnostics
    RuntimeAttr? runtime_attr_train_model
    RuntimeAttr? runtime_attr_shard_bed
    RuntimeAttr? runtime_attr_apply_mu_model
    RuntimeAttr? runtime_attr_merge_beds
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
  scatter ( contig in contigs ) {

    # Infer pairs BED file path
    String pairs_filename = select_first([pairs_bed_prefix, prefix + ".pairs.eigen"])
    File pairs_bed = pairs_bucket + "/" + pairs_filename + "." + contig + ".bed.gz" 
    File pairs_bed_idx = pairs_bed + ".tbi" 

    # Step 1. Intersect CNVs with bin-pairs
    call IntersectSVs {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        pairs_bed=pairs_bed,
        pairs_bed_idx=pairs_bed_idx,
        contig=contig,
        prefix="~{prefix}.~{cnv}",
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_override=runtime_attr_intersect_svs
    }

    # Step 2. Apply training mask
    call Utils.ApplyExclusionBED as ApplyTrainingMask {
      input:
        inbed=IntersectSVs.pairs_w_counts,
        exbed=training_mask,
        prefix=basename(IntersectSVs.pairs_w_counts, ".bed.gz") + ".training",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_apply_training_mask
    }
  }

  # Step 3. Train mutation rate model
  call TrainModel {
    input:
      training_beds=ApplyTrainingMask.filtered_bed,
      model=model,
      contigs_fai=contigs_fai,
      prefix="~{prefix}.~{cnv}",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_train_model
  }

  # Step 4. Apply mutation rate model to all bin-pairs, parallelized per chromosome
  # Prepare inputs for scatter
  Array[Pair[File, File]] predict_mu_beds_and_indexes = zip(IntersectSVs.pairs_w_counts, IntersectSVs.pairs_w_counts_idx)
  Array[Pair[String, Pair[File, File]]] predict_mu_inputs = zip(contigs, predict_mu_beds_and_indexes)

  # Scatter over chromosomes and apply mu model
  scatter ( inputs in predict_mu_inputs) {
    call ShardAndPredictMu.ShardAndPredictMu as PredictMu {
      input:
        bed=inputs.right.left,
        bed_idx=inputs.right.right,
        trained_model=TrainModel.trained_model,
        shard_size=shard_size_apply_mu,
        contig=inputs.left,
        prefix="~{prefix}.~{cnv}",
        athena_docker=athena_docker,
        runtime_attr_shard_bed=runtime_attr_shard_bed,
        runtime_attr_apply_mu_model=runtime_attr_apply_mu_model,
        runtime_attr_merge_beds=runtime_attr_merge_beds
    }
  }

  # Run diagnostics, if optioned
  if ( run_diagnostics ) {

    # Get stats of bin-pairs before and after applying training blacklist
    call GetPairDiagnostics {
      input:
        all_pairs=IntersectSVs.pairs_w_counts,
        training_pairs=ApplyTrainingMask.filtered_bed,
        cnv=cnv,
        prefix=prefix,
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_diagnostics
    }

    # TODO: plot training stats & model calibration from step 3

    # Tar all diagnostics for convenience
    call Utils.MakeTarball as MergeDiagnostics {
      input:
        files_to_tar=flatten([GetPairDiagnostics.all_outputs]),
        tarball_prefix="~{prefix}.~{cnv}.TrainMuModel.diagnostics",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_diagnostics
    }
  }

  output {
    Array[File] pairs_w_mu = PredictMu.pairs_w_mu
    Array[File] pairs_w_mu_idx = PredictMu.pairs_w_mu_idx
    File? diagnostics = MergeDiagnostics.tarball
  }
}


# Intersect SV breakpoints vs 2D bin-pairs for a single chromosome
task IntersectSVs {
  input {
    String vcf #VCF passed as string to allow for remote tabixing without localizing entire VCF
    File vcf_idx
    File pairs_bed
    File pairs_bed_idx
    String contig
    String prefix
    
    String athena_cloud_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 250,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Localize variants from contig
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{prefix}.~{contig}.svs.vcf.gz

    # Intersect variants and pairs with athena
    athena_cmd="athena count-sv --bin-format 2D --comparison breakpoint"
    athena_cmd="$athena_cmd --probabilities --bgzip ~{pairs_bed}"
    athena_cmd="$athena_cmd ~{prefix}.~{contig}.svs.vcf.gz"
    athena_cmd="$athena_cmd ~{prefix}.pairs.wCounts.~{contig}.bed.gz"
    echo -e "Now intersecting SVs and bins using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ~{prefix}.pairs.wCounts.~{contig}.bed.gz
  }

  output {
    File pairs_w_counts = "~{prefix}.pairs.wCounts.~{contig}.bed.gz"
    File pairs_w_counts_idx = "~{prefix}.pairs.wCounts.~{contig}.bed.gz.tbi"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: athena_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


# Collect diagnostic stats and plots for bin-pairs before and after filtering to training set
task GetPairDiagnostics {
  input {
    Array[File] all_pairs
    Array[File] training_pairs
    String cnv
    String prefix

    String dsmap_r_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 100,
    boot_disk_gb: 20,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Merge bin-pairs and keep first four columns only
    zcat ~{all_pairs[0]} | sed -n '1p' | cut -f1-4 > pairs.bed
    zcat ~{sep=" " all_pairs} | grep -ve '^#' | cut -f1-4 >> pairs.bed
    bgzip -f pairs.bed
    tabix -f pairs.bed.gz

    # Merge training bin-pairs and keep first four columns only
    zcat ~{training_pairs[0]} | sed -n '1p' | cut -f1-4 > training.bed
    zcat ~{sep=" " training_pairs} | grep -ve '^#' | cut -f1-4 >> training.bed
    bgzip -f training.bed
    tabix -f training.bed.gz

    # Get diagnostics
    mkdir outputs/
    /opt/dsmap/scripts/mu/get_pair_diagnostics.R \
      --cnv ~{cnv} \
      pairs.bed.gz \
      training.bed.gz \
      outputs/~{prefix}
  }

  output {
    Array[File] all_outputs = glob("outputs/~{prefix}*")
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: dsmap_r_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


# Train a CNV mutation rate model
task TrainModel {
  input {
    Array[File] training_beds
    String model
    File contigs_fai
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 8,
    disk_gb: 250,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Build list of training BEDs per contig
    while read contig; do
      find / -name "*.$contig.training.bed.gz" \
      | awk -v OFS="\t" -v contig="$contig" '{ print contig, $0 }'
    done < <( cut -f1 ~{contigs_fai} ) \
    > training_beds.tsv

    # Train model
    athena mu-train \
      --training-data training_beds.tsv \
      --model ~{model} \
      --model-outfile ~{prefix}.~{model}.trained.pkl \
      --stats-outfile ~{prefix}.~{model}.training_stats.tsv \
      --calibration-outfile ~{prefix}.~{model}.calibration.tsv
    gzip -f ~{prefix}.~{model}.calibration.tsv
  >>>

  output {
    File trained_model = "~{prefix}.~{model}.trained.pkl"
    File stats_tsv = "~{prefix}.~{model}.training_stats.tsv"
    File calibration_tsv = "~{prefix}.~{model}.calibration.tsv.gz"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: athena_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

