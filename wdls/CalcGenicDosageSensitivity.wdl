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

import "Utils.wdl"
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

    # Diagnostics options
    Boolean run_diagnostics = true

    # Dockers
    String athena_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_filter_gtf
    RuntimeAttr? runtime_attr_query_mu
    RuntimeAttr? runtime_attr_count_cnvs
    RuntimeAttr? runtime_attr_merge_data
    RuntimeAttr? runtime_attr_plot_mu_hist
    RuntimeAttr? runtime_attr_merge_diagnostics
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
  scatter ( contig in contigs ) {

    # Step 1. Filter GTF to coding sequences and gene bodies
    call FilterGtf {
      input:
        gtf=gtf,
        gtf_idx=gtf_idx,
        contig=contig,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_filter_gtf
    }

    File del_mu_bed = mu_bucket + "/" + mu_bed_prefix + ".DEL." + contig + ".mu.bed.gz"
    File del_mu_bed_idx = del_mu_bed + ".tbi"
    File dup_mu_bed = mu_bucket + "/" + mu_bed_prefix + ".DUP." + contig + ".mu.bed.gz"
    File dup_mu_bed_idx = dup_mu_bed + ".tbi" 

    # Step 2a. Compute mutation rates for all CDS-overlapping deletions per gene
    call QueryMuForGtf as QueryMuCodingDel {
      input:
        gtf=FilterGtf.coding_gtf,
        gtf_idx=FilterGtf.coding_gtf_idx,
        mu_bed=del_mu_bed,
        mu_bed_idx=del_mu_bed_idx,
        athena_query_options=[""],
        prefix=basename(FilterGtf.coding_gtf, ".gtf.gz") + ".DEL.CDS",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 2b. Count CDS-overlapping deletions per gene
    call CountCnvs as CountCodingDel {
      input:
        vcf=del_vcf,
        vcf_idx=del_vcf_idx,
        gtf=FilterGtf.coding_gtf,
        gtf_idx=FilterGtf.coding_gtf_idx,
        athena_countsv_options=[""],
        prefix=basename(FilterGtf.coding_gtf, ".gtf.gz") + ".DEL.CDS.counts",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }

    # Step 3a. Compute mutation rates for all CDS-overlapping duplications per gene
    call QueryMuForGtf as QueryMuCodingDup {
      input:
        gtf=FilterGtf.coding_gtf,
        gtf_idx=FilterGtf.coding_gtf_idx,
        mu_bed=dup_mu_bed,
        mu_bed_idx=dup_mu_bed_idx,
        athena_query_options=[""],
        prefix=basename(FilterGtf.coding_gtf, ".gtf.gz") + ".DUP.CDS",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 3b. Count CDS-overlapping duplications per gene
    call CountCnvs as CountCodingDup {
      input:
        vcf=dup_vcf,
        vcf_idx=dup_vcf_idx,
        gtf=FilterGtf.coding_gtf,
        gtf_idx=FilterGtf.coding_gtf_idx,
        athena_countsv_options=[""],
        prefix=basename(FilterGtf.coding_gtf, ".gtf.gz") + ".DUP.CDS.counts",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }

    # Step 4a. Compute mutation rates for all copy-gain duplications per gene
    call QueryMuForGtf as QueryMuCopyGainDup {
      input:
        gtf=FilterGtf.genes_gtf,
        gtf_idx=FilterGtf.genes_gtf_idx,
        mu_bed=dup_mu_bed,
        mu_bed_idx=dup_mu_bed_idx,
        athena_query_options=["--fraction 1.0"],
        prefix=basename(FilterGtf.genes_gtf, ".gtf.gz") + ".DUP.CG",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 4b. Count copy-gain duplications per gene
    call CountCnvs as CountCopyGainDup {
      input:
        vcf=dup_vcf,
        vcf_idx=dup_vcf_idx,
        gtf=FilterGtf.genes_gtf,
        gtf_idx=FilterGtf.genes_gtf_idx,
        athena_countsv_options=["--fraction 1.0"],
        prefix=basename(FilterGtf.genes_gtf, ".gtf.gz") + ".DUP.CG.counts",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }
  }

  # Step 2c. Merge and analyze outputs from 2a & 2b
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeCodingDelData {
    input:
      mu_tsvs=QueryMuCodingDel.mu_tsv,
      counts_tsvs=CountCodingDel.counts_tsv,
      prefix=prefix + ".DEL.CDS",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }
  
  # Step 3c. Merge and analyze outputs from 3a & 3b
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeCodingDupData {
    input:
      mu_tsvs=QueryMuCodingDup.mu_tsv,
      counts_tsvs=CountCodingDup.counts_tsv,
      prefix=prefix + ".DUP.CDS",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }

  # Step 4c. Merge and analyze outputs from 4a & 4b
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeCopyGainDupData {
    input:
      mu_tsvs=QueryMuCopyGainDup.mu_tsv,
      counts_tsvs=CountCopyGainDup.counts_tsv,
      prefix=prefix + ".DUP.CG",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }

  # Gather diagnostics, if optioned
  if ( run_diagnostics ) {

    # Plot mutation rates per gene as a sanity check
    call Utils.PlotMuHist as PlotCodingDelMu {
      input:
        mu_tsv=MergeCodingDelData.merged_tsv,
        cnv="DEL",
        x_title="Coding deletions per allele per generation",
        y_title="Genes",
        out_prefix=prefix + ".DEL.CDS",
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuHist as PlotCodingDupMu {
      input:
        mu_tsv=MergeCodingDupData.merged_tsv,
        cnv="DUP",
        x_title="Coding duplications per allele per generation",
        y_title="Genes",
        out_prefix=prefix + ".DUP.CDS",
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuHist as PlotCopyGainDupMu {
      input:
        mu_tsv=MergeCopyGainDupData.merged_tsv,
        cnv="DUP",
        x_title="Whole-gene duplications per allele per generation",
        y_title="Genes",
        out_prefix=prefix + ".DUP.CG",
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }

    # Tar diagnostics, for convenience
    call Utils.MakeTarball as MergeDiagnostics {
      input:
        files_to_tar=[PlotCodingDelMu.mu_hist, PlotCodingDupMu.mu_hist, PlotCopyGainDupMu.mu_hist],
        tarball_prefix="~{prefix}.CalcGenicDosageSensitivity.diagnostics",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_merge_diagnostics
    }
  }

  output {
    File coding_del_data_tsv = MergeCodingDelData.merged_tsv
    File coding_dup_data_tsv = MergeCodingDupData.merged_tsv
    File copy_gain_dup_data_tsv = MergeCopyGainDupData.merged_tsv
    File? diagnostics = MergeDiagnostics.tarball
  }
}


# Filter GTF to a single chromosome and subsets of CDSs or gene bodies
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
    set -euo pipefail

    # Subset GTF to coding sequences
    tabix -h ~{gtf} ~{contig} | awk '{ if ($3 == "CDS") print $0 }' | bgzip -c \
    > ~{gtf_prefix}.coding.~{contig}.gtf.gz
    tabix -f ~{gtf_prefix}.coding.~{contig}.gtf.gz

    # Subset GTF to gene bodies
    tabix -h ~{gtf} ~{contig} | awk '{ if ($3 == "gene") print $0 }' | bgzip -c \
    > ~{gtf_prefix}.genes.~{contig}.gtf.gz
    tabix -f ~{gtf_prefix}.genes.~{contig}.gtf.gz
  >>>

  output {
    File coding_gtf = "~{gtf_prefix}.coding.~{contig}.gtf.gz"
    File coding_gtf_idx = "~{gtf_prefix}.coding.~{contig}.gtf.gz.tbi"
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
    set -euo pipefail

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
    String prefix

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
    set -euo pipefail

    # Count SVs
    athena_cmd="athena count-sv --query-format gtf --outfile ~{prefix}.tsv.gz"
    athena_cmd="$athena_cmd --bgzip  ~{sep=' ' athena_countsv_options}"
    athena_cmd="$athena_cmd ~{vcf} ~{gtf}"
    echo -e "Now counting SVs using command:\n$athena_cmd"
    eval $athena_cmd
  >>>

  output {
    File counts_tsv = "~{prefix}.tsv.gz"
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


# Merge gene-level mutation rates and CNV counts
# TODO: add analysis component to this. Currently just merges outputs
task MergeMuAndCounts {
  input {
    Array[File] mu_tsvs
    Array[File] counts_tsvs
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override    
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(10 * size(flatten([mu_tsvs, counts_tsvs]), "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Merge mutation rates
    zcat ~{mu_tsvs[0]} | sed -n '1p' > mu.header
    zcat ~{sep=" " mu_tsvs} | grep -ve '^#' | sort -k1,1 | cat mu.header - | gzip -c \
    > ~{prefix}.mu.tsv.gz

    # Merge counts
    zcat ~{counts_tsvs[0]} | sed -n '1p' > counts.header
    zcat ~{sep=" " counts_tsvs} | grep -ve '^#' | sort -k1,1 | cat counts.header - | gzip -c \
    > ~{prefix}.counts.tsv.gz

    # Join mutation rates and counts
    join -j 1 -t $'\t' \
      <( zcat ~{prefix}.mu.tsv.gz ) \
      <( zcat ~{prefix}.counts.tsv.gz ) \
    | gzip -c \
    > ~{prefix}.mu_and_counts.tsv.gz
  >>>

  output {
    File merged_tsv = "~{prefix}.mu_and_counts.tsv.gz"
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
