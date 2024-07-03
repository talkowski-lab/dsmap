#######################
#    DSMap Project    #
#######################
#
# overlapSV.wdl
# only overlapSV task from TrainMuModel

version 1.0

import "Utils.wdl"
import "Structs.wdl"

workflow TrainMuModel {
  input {
    # General inputs
    String vcf
    File vcf_idx
    String pairs_bucket
    String? pairs_bed_prefix
    File contigs_fai
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
  }
  output {
    Array[File] pairs_w_counts = IntersectSVs.pairs_w_counts
    Array[File] pairs_w_counts_idx = IntersectSVs.pairs_w_counts_idx
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

  RuntimeAttr default_attr = {
    "cpu_cores": 1, 
    "mem_gb": 64,
    "disk_gb": 50,
    "boot_disk_gb": 10,
    "preemptible_tries": 1,
    "max_retries": 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Localize variants from contig
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{prefix}.~{contig}.svs.vcf.gz

    # Intersect variants and pairs with athena
    athena_cmd="athena count-sv --query-format pairs --comparison breakpoint"
    athena_cmd="$athena_cmd --probabilities --bgzip"
    athena_cmd="$athena_cmd --outfile ~{prefix}.pairs.wCounts.~{contig}.bed.gz"
    athena_cmd="$athena_cmd ~{prefix}.~{contig}.svs.vcf.gz ~{pairs_bed}"
    echo -e "Now intersecting SVs and bins using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ~{prefix}.pairs.wCounts.~{contig}.bed.gz

    # Count pairs (used later in TrainMu)
    zcat ~{pairs_bed} | grep -ve '^#' | cut -f1 | wc -l > n_pairs.txt
  }

  output {
    File pairs_w_counts = "~{prefix}.pairs.wCounts.~{contig}.bed.gz"
    File pairs_w_counts_idx = "~{prefix}.pairs.wCounts.~{contig}.bed.gz.tbi"
    Int n_pairs = read_int("n_pairs.txt")
  }
  
  runtime {
    cpu: runtime_attr["cpu_cores"]
    memory: runtime_attr["mem_gb"] + " GiB"
    disks: "local-disk " + runtime_attr["disk_gb"] + " HDD"
    bootDiskSizeGb: runtime_attr["boot_disk_gb"]
    preemptible: runtime_attr["preemptible_tries"]
    maxRetries: runtime_attr["max_retries"]
    docker: athena_cloud_docker
  }
}
