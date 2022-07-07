#######################
#    DSMap Project    #
#######################
#
# BinAndAnnotateGenome.wdl
#
# Bin and annotate a genome or set of chromosomes
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"
import "AnnotateBinsSingleChrom.wdl"
import "MakeAndAnnotatePairsSingleChrom.wdl"
import "DecompAnnotationsSingleChrom.wdl"


workflow BinAndAnnotateGenome {
  input {
    # General inputs
    File contigs_fai
    String prefix
    String? ref_build = "hg38"
    File? ref_fasta

    # Bin inputs
    File bin_exclusion_mask
    Int bin_size
    Int bins_per_shard
    File? bin_annotations_list_localize
    File? bin_annotations_list_remote
    File? bin_annotations_list_ucsc
    File? snv_mutrates_tsv
    
    # Pair inputs
    File? pair_exclusion_mask
    Int max_pair_distance
    Int bins_per_pair_shard
    File? pair_annotations_list_localize
    File? pair_annotations_list_remote
    File? pair_annotations_list_ucsc

    # PCA inputs
    Boolean decompose_features = true
    Boolean visualize_features_before_pca = false
    File? feature_transformations_tsv
    Int? pairs_for_pca
    Float? pca_min_variance
    Int? max_pcs
    Int? pairs_per_shard_apply_pca

    # Dockers
    String athena_docker
    String athena_cloud_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_make_bins
    RuntimeAttr? runtime_attr_chrom_shard
    RuntimeAttr? runtime_attr_annotate_bins
    RuntimeAttr? runtime_attr_merge_annotated_bins
    RuntimeAttr? runtime_attr_make_pairs
    RuntimeAttr? runtime_attr_annotate_pairs
    RuntimeAttr? runtime_attr_sample_pairs
    RuntimeAttr? runtime_attr_merge_pairs
    RuntimeAttr? runtime_attr_visualize_features
    RuntimeAttr? runtime_attr_learn_pca
    RuntimeAttr? runtime_attr_apply_pca
  }

<<<<<<< HEAD
=======

>>>>>>> main
  # Step 1. Create all 1D bins
  call MakeBins {
    input:
      contigs_fai=contigs_fai,
      bin_size=bin_size,
      bin_exclusion_mask=bin_exclusion_mask,
      prefix=prefix,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_make_bins
  }

  # Prior to parallelizing per chromosome, must determine number of pairs to sample
  if ( decompose_features ) { 
    call CalcPairsPerChrom {
      input:
        contigs_fai=contigs_fai,
        pairs_for_pca=pairs_for_pca,
        prefix=prefix,
        athena_docker=athena_docker
    }
  }

  Array[Array[String]] contigs = read_tsv(select_first([CalcPairsPerChrom.updated_fai, contigs_fai]))

  # Process each chromosome in parallel
  scatter ( contig in contigs ) {

    # Step 2. Annotate 1D bins per contig
    call AnnotateBinsSingleChrom.AnnotateBinsSingleChrom as AnnotateBins {
      input:
        bins=MakeBins.bins_bed,
        bins_idx=MakeBins.bins_bed_idx,
        bedtools_genome_file=MakeBins.bedtools_genome_file,
        contig=contig[0],
        shard_size=bins_per_shard,
        prefix=prefix,
        bin_annotations_list_localize=bin_annotations_list_localize,
        bin_annotations_list_remote=bin_annotations_list_remote,
        bin_annotations_list_ucsc=bin_annotations_list_ucsc,
        ref_build=ref_build,
        ref_fasta=ref_fasta,
        snv_mutrates_tsv=snv_mutrates_tsv,
        athena_docker=athena_docker,
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_chrom_shard=runtime_attr_chrom_shard,
        runtime_attr_annotate_bins=runtime_attr_annotate_bins,
        runtime_attr_merge_annotated_bins=runtime_attr_merge_annotated_bins
    }

    # Step 3. Pair 2D bins and add 2D annotations
    call MakeAndAnnotatePairsSingleChrom.MakeAndAnnotatePairsSingleChrom as MakeAndAnnotatePairs {
      input:
        bins=AnnotateBins.annotated_bins,
        bins_idx=AnnotateBins.annotated_bins_idx,
        bedtools_genome_file=MakeBins.bedtools_genome_file,
        pair_exclusion_mask=pair_exclusion_mask,
        contig=contig[0],
        max_pair_distance=max_pair_distance,
        shard_size=bins_per_pair_shard,
        prefix=prefix,
        pair_annotations_list_localize=pair_annotations_list_localize,
        pair_annotations_list_remote=pair_annotations_list_remote,
        pair_annotations_list_ucsc=pair_annotations_list_ucsc,
        ref_build=ref_build,
        ref_fasta=ref_fasta,
        bin_size=bin_size,
        sample_pairs_for_pca=decompose_features,
        pairs_to_sample_for_pca=contig[2],
        athena_docker=athena_docker,
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_chrom_shard=runtime_attr_chrom_shard,
        runtime_attr_make_pairs=runtime_attr_make_pairs,
        runtime_attr_annotate_pairs=runtime_attr_annotate_pairs,
        runtime_attr_sample_pairs=runtime_attr_sample_pairs,
        runtime_attr_merge_annotated_pairs=runtime_attr_merge_pairs
    }
  }

  # Steps 4-5 are only necessary to run if decompose_features is specified
  if ( decompose_features ) { 
<<<<<<< HEAD

    # Step 4.1. Merge subsampled pairs for PCA
    call Utils.MergeBEDs as MergePCAPairs {
      input:
        beds=select_all(MakeAndAnnotatePairs.downsampled_pairs),
        prefix="~{prefix}.annotated_pairs.downsampled",
        beds_are_bgzipped=true,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_merge_pairs
    }

    # [Optional] Step 4.2. Visualize distributions of all features prior to PCA
    if ( visualize_features_before_pca ) {
      call Utils.VisualizeFeatures as VisualizeFeatures {
        input:
          bed=MergePCAPairs.merged_bed,
          transformations_tsv=feature_transformations_tsv,
          prefix=prefix,
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_visualize_features
      }
    }

    # Step 4.3. Learn transformation
    call LearnPCA {
      input:
        sampled_pairs=MergePCAPairs.merged_bed,
        sampled_pairs_idx=MergePCAPairs.merged_bed_idx,
        pca_min_variance=pca_min_variance,
        max_pcs=max_pcs,
        transformations_tsv=feature_transformations_tsv,
        prefix=prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_learn_pca
    }

    # Prior to Step 5, need to compose grouped inputs per chromosome
    Array[Pair[File, File]] pairs_with_idxs = zip(MakeAndAnnotatePairs.annotated_pairs, MakeAndAnnotatePairs.annotated_pairs_idx)
    Array[Pair[Pair[File, File], String]] inputs_for_decomp = zip(pairs_with_idxs, transpose(contigs)[0])

    # Step 5. Apply PCA transformation to 2D pairs (in parallel)
    scatter ( contig_inputs in inputs_for_decomp ) {

      # For readability, extract inputs here
      File chrom_pairs_bed = contig_inputs.left.left
      File chrom_pairs_bed_idx = contig_inputs.left.right
      String contig = contig_inputs.right

      # Apply PCA transformation
      call DecompAnnotationsSingleChrom.DecompAnnotationsSingleChrom as DecompAnnosPerChrom {
        input:
          pairs=chrom_pairs_bed,
          pairs_idx=chrom_pairs_bed_idx,
          pca_model=LearnPCA.pca_model,
          contig=contig,
          shard_size=pairs_per_shard_apply_pca,
          prefix="~{prefix}.pairs",
          athena_docker=athena_docker,
          runtime_attr_chrom_shard=runtime_attr_chrom_shard,
          runtime_attr_apply_pca=runtime_attr_apply_pca,
          runtime_attr_merge_pairs=runtime_attr_merge_pairs
      }
    }
  }

  output {

    Array[File] annotated_pairs = MakeAndAnnotatePairs.annotated_pairs
    Array[File] annotated_pairs_idx = MakeAndAnnotatePairs.annotated_pairs_idx

    Array[File]? decomped_pairs = DecompAnnosPerChrom.decomped_bed
    Array[File]? decomped_pairs_idx = DecompAnnosPerChrom.decomped_bed_idx

    File? feature_distribs_pre_pca = VisualizeFeatures.feature_hists

=======

    # Step 4.1. Merge subsampled pairs for PCA
    call Utils.MergeBEDs as MergePCAPairs {
      input:
        beds=select_all(MakeAndAnnotatePairs.downsampled_pairs),
        prefix="~{prefix}.annotated_pairs.downsampled",
        beds_are_bgzipped=true,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_merge_pairs
    }

    # [Optional] Step 4.2. Visualize distributions of all features prior to PCA
    if ( visualize_features_before_pca ) {
      call Utils.VisualizeFeatures as VisualizeFeatures {
        input:
          bed=MergePCAPairs.merged_bed,
          transformations_tsv=feature_transformations_tsv,
          prefix=prefix,
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_visualize_features
      }
    }

    # Step 4.3. Learn transformation
    call LearnPCA {
      input:
        sampled_pairs=MergePCAPairs.merged_bed,
        sampled_pairs_idx=MergePCAPairs.merged_bed_idx,
        pca_min_variance=pca_min_variance,
        max_pcs=max_pcs,
        transformations_tsv=feature_transformations_tsv,
        prefix=prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_learn_pca
    }

    # Prior to Step 5, need to compose grouped inputs per chromosome
    Array[Pair[File, File]] pairs_with_idxs = zip(MakeAndAnnotatePairs.annotated_pairs, MakeAndAnnotatePairs.annotated_pairs_idx)
    Array[Pair[Pair[File, File], String]] inputs_for_decomp = zip(pairs_with_idxs, transpose(contigs)[0])

    # Step 5. Apply PCA transformation to 2D pairs (in parallel)
    scatter ( contig_inputs in inputs_for_decomp ) {

      # For readability, extract inputs here
      File chrom_pairs_bed = contig_inputs.left.left
      File chrom_pairs_bed_idx = contig_inputs.left.right
      String contig = contig_inputs.right

      # Apply PCA transformation
      call DecompAnnotationsSingleChrom.DecompAnnotationsSingleChrom as DecompAnnosPerChrom {
        input:
          pairs=chrom_pairs_bed,
          pairs_idx=chrom_pairs_bed_idx,
          pca_model=LearnPCA.pca_model,
          contig=contig,
          shard_size=pairs_per_shard_apply_pca,
          prefix="~{prefix}.pairs",
          athena_docker=athena_docker,
          runtime_attr_chrom_shard=runtime_attr_chrom_shard,
          runtime_attr_apply_pca=runtime_attr_apply_pca,
          runtime_attr_merge_pairs=runtime_attr_merge_pairs
      }
    }
  }

  output {

    Array[File] annotated_pairs = MakeAndAnnotatePairs.annotated_pairs
    Array[File] annotated_pairs_idx = MakeAndAnnotatePairs.annotated_pairs_idx

    Array[File]? decomped_pairs = DecompAnnosPerChrom.decomped_bed
    Array[File]? decomped_pairs_idx = DecompAnnosPerChrom.decomped_bed_idx

    File? feature_distribs_pre_pca = VisualizeFeatures.feature_hists

>>>>>>> main
    File? eigenfeature_stats = LearnPCA.pc_stats

  }
}


# Create sequential uniform bins for all chromosomes
task MakeBins {
  input {
    File contigs_fai
    Int bin_size
    File bin_exclusion_mask
    String prefix
    
    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 25,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Create .genome file from contigs_fai
    cut -f1-2 ~{contigs_fai} > contigs.genome

    # Make bins & tabix output BED
    athena make-bins \
      --exclusion-list-all ~{bin_exclusion_mask} \
      --buffer ~{bin_size} \
      --include-chroms $( cut -f1 contigs.genome | paste -s -d, ) \
      --bgzip \
      contigs.genome \
      ~{bin_size} \
      ~{prefix}.bins.bed.gz
    tabix -f ~{prefix}.bins.bed.gz
  }

  output {
    File bins_bed = "~{prefix}.bins.bed.gz"
    File bins_bed_idx = "~{prefix}.bins.bed.gz.tbi"
    File bedtools_genome_file = "contigs.genome"
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


# Evaluate number of pairs to sample per chromosome corresponding to length
task CalcPairsPerChrom {
  input {
    File contigs_fai
    Int pairs_for_pca = 100000
    String prefix

    String athena_docker
  }
  
  runtime {
    cpu: 1
    memory: "3 GiB"
    disks: "local-disk 50 HDD"
    bootDiskSizeGb: 10
    docker: athena_docker
    preemptible: 3
    maxRetries: 1
  }

  command <<<
<<<<<<< HEAD
    set -eu -o pipefail

    total_bp=$( awk -v FS="\t" '{ sum+=$2 }END{ printf "%.0f", sum }' ~{contigs_fai} )
=======
    set -euo pipefail

    total_bp=$( awk -v FS="\t" '{ sum+=$2 }END{ print sum }' ~{contigs_fai} )
>>>>>>> main
    bp_per_pair=$(( ( $total_bp + ~{pairs_for_pca} - 1 ) / ~{pairs_for_pca} ))

    # Produces a revised .fai file with chromosome, length, and number of pairs to sample
    awk -v bp_per_pair="$bp_per_pair" -v FS="\t" \
      '{ printf "%s\t%s\t%0.0f\n", $1, $2, ($2 + bp_per_pair - 1) / bp_per_pair }' \
      ~{contigs_fai} \
    > ~{prefix}.updated.fai
  >>>

  output {
    File updated_fai = "~{prefix}.updated.fai"
  }
}


# Learn PCA transformation to decompose features for all pairs
task LearnPCA {
  input {
    File sampled_pairs
    File sampled_pairs_idx
    String prefix

    File? transformations_tsv
    Float? pca_min_variance
    Int? max_pcs

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
<<<<<<< HEAD
    cpu_cores: 2, 
    mem_gb: 4,
    disk_gb: 10 + (10 * ceil(size(sampled_pairs, "GB"))),
=======
    cpu_cores: 4, 
    mem_gb: 16,
    disk_gb: 250,
>>>>>>> main
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Create option for feature transformation prior to model learning
    athena_options=""
    if [ "~{defined(transformations_tsv)}" == "true" ]; then
      athena_options="$athena_options --transformations-tsv ~{transformations_tsv}"
    fi
    if [ "~{defined(pca_min_variance)}" == "true" ]; then
      athena_options="$athena_options --min-variance ~{pca_min_variance}"
    fi
    if [ "~{defined(max_pcs)}" == "true" ]; then
      athena_options="$athena_options --max-components ~{max_pcs}"
    fi

    # Learn PCA transformation with athena
    athena_cmd="athena eigen-bins $athena_options --bgzip"
<<<<<<< HEAD
    athena_cmd="$athena_cmd --parameters-outfile ~{prefix}.pca_model.pkl"
=======
    athena_cmd="$athena_cmd --parameters-outfile ~{prefix}.pca_model.pickle"
>>>>>>> main
    athena_cmd="$athena_cmd --whiten --fill-missing mean"
    athena_cmd="$athena_cmd --stats ~{prefix}.pca_stats.txt ~{sampled_pairs}"
    echo -e "Now learning PCA decomposition using command:\n$athena_cmd"
    eval $athena_cmd
  }

  output {
<<<<<<< HEAD
    File pca_model = "~{prefix}.pca_model.pkl"
=======
    File pca_model = "~{prefix}.pca_model.pickle"
>>>>>>> main
    File pc_stats = "~{prefix}.pca_stats.txt"
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


