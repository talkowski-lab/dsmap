#######################
#    DSMap Project    #
#######################
#
# CalcGenicDosageSensitivity.wdl
#
# Compute dosage sensitivity statistics for all genes in a single SV dataset
#
# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Structs.wdl"


workflow CalcGenicDosageSensitivity {
  input {
    # General inputs
    File del_vcf
    File del_vcf_idx
    File dup_vcf
    File dup_vcf_idx
    File gtf
    File gtf_idx
    String mu_bucket
    String mu_bed_prefix
    File contigs_fai
    String prefix

    # Dockers
    String athena_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_filter_gtf
    RuntimeAttr? runtime_attr_query_mu
    RuntimeAttr? runtime_attr_count_cnvs
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
  scatter ( contig in contigs ) {

    # Step 1. Filter GTF to exons and gene bodies
    call FilterGtf {
      input:
        gtf=gtf,
        gtf_idx=gtf_idx,
        contig=contig,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_filter_gtf
    }

    # Step 2a. Compute mutation rates for all exon-overlapping deletions per gene
    File del_mu_bed = mu_bucket + "/" + mu_bed_prefix + ".DEL." + contig + ".mu.bed.gz"
    File del_mu_bed_idx = del_mu_bed + ".tbi"
    call QueryMuForGtf as QueryMuHaplo {
      input:
        gtf=FilterGtf.exons_gtf,
        gtf_idx=FilterGtf.exons_gtf_idx,
        mu_bed=del_mu_bed,
        mu_bed_idx=del_mu_bed_idx,
        athena_query_options=[""],
        prefix=basename(FilterGtf.exons_gtf, ".gtf.gz") + ".DEL",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 2b. Count exon-overlapping deletions per gene
    call CountCnvs as CountDel {
      input:
        vcf=del_vcf,
        vcf_idx=del_vcf_idx,
        gtf=FilterGtf.exons_gtf,
        gtf_idx=FilterGtf.exons_gtf_idx,
        athena_countsv_options=[""],
        output_prefix=basename(FilterGtf.exons_gtf, ".gtf.gz") + ".DEL.counts",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }

    # Step 3a. Compute mutation rates for all copy-gain duplications per gene
    File dup_mu_bed = mu_bucket + "/" + mu_bed_prefix + ".DUP." + contig + ".mu.bed.gz"
    File dup_mu_bed_idx = dup_mu_bed + ".tbi" 
    call QueryMuForGtf as QueryMuTriplo {
      input:
        gtf=FilterGtf.genes_gtf,
        gtf_idx=FilterGtf.genes_gtf_idx,
        mu_bed=dup_mu_bed,
        mu_bed_idx=dup_mu_bed_idx,
        athena_query_options=["--fraction 1.0"],
        prefix=basename(FilterGtf.genes_gtf, ".gtf.gz") + ".DUP",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 3b. Count copy-gain duplications per gene
    call CountCnvs as CountDup {
      input:
        vcf=dup_vcf,
        vcf_idx=dup_vcf_idx,
        gtf=FilterGtf.genes_gtf,
        gtf_idx=FilterGtf.genes_gtf_idx,
        athena_countsv_options=["--fraction 1.0"],
        output_prefix=basename(FilterGtf.genes_gtf, ".gtf.gz") + ".DUP.counts",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }
  }

  # Step 2c. Merge and analyze outputs from 2a & 2b
  # TODO: implement this

  # Step 3c. Merge and analyze outputs from 3a & 3b
  # TODO: implement this

  output {}
}


# Filter GTF to a single chromosome and subsets of exons or gene bodies
task FilterGtf {
  input {
    File gtf
    File gtf_idx
    String contig

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  String gtf_prefix = basename(gtf, ".gtf.gz")

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size([gtf], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -eu -o pipefail

    # Subset GTF to exons
    tabix -h ~{gtf} ~{contig} | awk '{ if ($3 == "exon") print $0 }' | bgzip -c \
    > ~{gtf_prefix}.exons.~{contig}.gtf.gz
    tabix -f ~{gtf_prefix}.exons.~{contig}.gtf.gz

    # Subset gtf to gene bodies
    tabix -h ~{gtf} ~{contig} | awk '{ if ($3 == "gene") print $0 }' | bgzip -c \
    > ~{gtf_prefix}.genes.~{contig}.gtf.gz
    tabix -f ~{gtf_prefix}.genes.~{contig}.gtf.gz
  >>>

  output {
    File exons_gtf = "~{gtf_prefix}.exons.~{contig}.gtf.gz"
    File exons_gtf_idx = "~{gtf_prefix}.exons.~{contig}.gtf.gz.tbi"
    File genes_gtf = "~{gtf_prefix}.genes.~{contig}.gtf.gz"
    File genes_gtf_idx = "~{gtf_prefix}.genes.~{contig}.gtf.gz.tbi"
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


# Compute mutation rates for an input GTF
task QueryMuForGtf {
  input {
    File gtf
    File? gtf_idx
    File mu_bed
    File mu_bed_idx
    Array[String]? athena_query_options = [""]
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 10 + ceil(2 * size([gtf, mu_bed], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -eu -o pipefail

    # Query mutation rates
    athena_cmd="athena mu-query --gzip -o ~{prefix}.mu.tsv.gz"
    athena_cmd="$athena_cmd ~{sep=' ' athena_query_options}"
    athena_cmd="$athena_cmd ~{mu_bed} ~{gtf}"
    echo -e "Now querying mutation rates using command:\n$athena_cmd"
    eval $athena_cmd
  }

  output {
    File mu_tsv = "~{prefix}.mu.tsv.gz"
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


# Count number of qualifying CNVs per gene
task CountCnvs {
  input {
    File vcf
    File vcf_idx
    File gtf
    File gtf_idx
    Array[String]? athena_countsv_options = [""]
    String output_prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size([gtf, vcf], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -eu -o pipefail

    # Count SVs
    athena_cmd="athena count-sv --query-format gtf --outfile ~{output_prefix}.tsv.gz"
    athena_cmd="$athena_cmd --bgzip  ~{sep=' ' athena_countsv_options}"
    athena_cmd="$athena_cmd ~{vcf} ~{gtf}"
    echo -e "Now counting SVs using command:\n$athena_cmd"
    eval $athena_cmd
  >>>

  output {
    File counts_tsv = "~{output_prefix}.tsv.gz"
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
