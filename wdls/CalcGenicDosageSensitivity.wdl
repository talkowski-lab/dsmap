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
# Contact: Lily Wang <lily_wang@hms.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"
import "CalcSegmentDosageSensitivity.wdl"


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
    String athena_cloud_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_filter_gtf
    RuntimeAttr? runtime_attr_filter_gtf_chrom
    RuntimeAttr? runtime_attr_query_mu
    RuntimeAttr? runtime_attr_expand_query
    RuntimeAttr? runtime_attr_count_cnvs
    RuntimeAttr? runtime_attr_format_counts
    RuntimeAttr? runtime_attr_merge_data
    RuntimeAttr? runtime_attr_plot_mu_hist
    RuntimeAttr? runtime_attr_merge_diagnostics
  }

  # Extract CDSs and gene bodies from GTF
  call FilterGtf {
    input:
      gtf=gtf,
      gtf_idx=gtf_idx,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_filter_gtf
  }

  # Calculate dosage sensitivity for whole gene copy-losses and copy-gains
  call CalcSegmentDosageSensitivity.CalcSegmentDosageSensitivity as CalcGeneCopyDosageSensitivity {
    input:
      del_vcf=del_vcf,
      del_vcf_idx=del_vcf_idx,
      dup_vcf=dup_vcf,
      dup_vcf_idx=dup_vcf_idx,
      query=FilterGtf.genes_gtf,
      query_idx=FilterGtf.genes_gtf_idx,
      mu_bucket=mu_bucket,
      mu_bed_prefix=mu_bed_prefix,
      contigs_fai=contigs_fai,
      prefix=prefix + ".gene",
      full_segment_overlap=true,
      run_diagnostics=run_diagnostics,
      athena_docker=athena_docker,
      athena_cloud_docker=athena_cloud_docker,
      dsmap_r_docker=dsmap_r_docker,
      runtime_attr_filter_query_chrom=runtime_attr_filter_gtf_chrom,
      runtime_attr_query_mu=runtime_attr_query_mu,
      runtime_attr_expand_query=runtime_attr_expand_query,
      runtime_attr_count_cnvs=runtime_attr_count_cnvs,
      runtime_attr_format_counts=runtime_attr_format_counts,
      runtime_attr_merge_data=runtime_attr_merge_data,
      runtime_attr_plot_mu_hist=runtime_attr_plot_mu_hist,
      runtime_attr_merge_diagnostics=runtime_attr_merge_diagnostics
  }

  # Calculate dosage sensitivity for whole gene copy-losses and copy-gains
  call CalcSegmentDosageSensitivity.CalcSegmentDosageSensitivity as CalcCodingDosageSensitivity {
    input:
      del_vcf=del_vcf,
      del_vcf_idx=del_vcf_idx,
      dup_vcf=dup_vcf,
      dup_vcf_idx=dup_vcf_idx,
      query=FilterGtf.coding_gtf,
      query_idx=FilterGtf.coding_gtf_idx,
      mu_bucket=mu_bucket,
      mu_bed_prefix=mu_bed_prefix,
      contigs_fai=contigs_fai,
      prefix=prefix + ".CDS",
      full_segment_overlap=false,
      run_diagnostics=run_diagnostics,
      athena_docker=athena_docker,
      athena_cloud_docker=athena_cloud_docker,
      dsmap_r_docker=dsmap_r_docker,
      runtime_attr_filter_query_chrom=runtime_attr_filter_gtf_chrom,
      runtime_attr_query_mu=runtime_attr_query_mu,
      runtime_attr_expand_query=runtime_attr_expand_query,
      runtime_attr_count_cnvs=runtime_attr_count_cnvs,
      runtime_attr_format_counts=runtime_attr_format_counts,
      runtime_attr_merge_data=runtime_attr_merge_data,
      runtime_attr_plot_mu_hist=runtime_attr_plot_mu_hist,
      runtime_attr_merge_diagnostics=runtime_attr_merge_diagnostics
  }

  # Gather diagnostics, if optioned
  if ( run_diagnostics ) {

    # Plot mutation rates per gene as a sanity check
    call Utils.PlotMuHist as PlotCodingDelMu {
      input:
        mu_tsv=CalcCodingDosageSensitivity.del_data_tsv,
        cnv="DEL",
        x_title="Coding deletions per allele per generation",
        y_title="Genes",
        out_prefix=prefix + ".DEL.CDS",
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuHist as PlotCodingDupMu {
      input:
        mu_tsv=CalcCodingDosageSensitivity.dup_data_tsv,
        cnv="DUP",
        x_title="Coding duplications per allele per generation",
        y_title="Genes",
        out_prefix=prefix + ".DUP.CDS",
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuHist as PlotCopyLossDelMu {
      input:
        mu_tsv=CalcGeneCopyDosageSensitivity.del_data_tsv,
        cnv="DEL",
        x_title="Whole-gene losses per allele per generation",
        y_title="Genes",
        out_prefix=prefix + ".DEL.CL",
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuHist as PlotCopyGainDupMu {
      input:
        mu_tsv=CalcGeneCopyDosageSensitivity.dup_data_tsv,
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
        files_to_tar=[PlotCodingDelMu.mu_hist, PlotCodingDupMu.mu_hist,
                      PlotCopyLossDelMu.mu_hist, PlotCopyGainDupMu.mu_hist],
        tarball_prefix="~{prefix}.CalcGenicDosageSensitivity.diagnostics",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_merge_diagnostics
    }
  }

  output {
    File coding_del_data_tsv = CalcCodingDosageSensitivity.del_data_tsv
    File coding_dup_data_tsv = CalcCodingDosageSensitivity.dup_data_tsv
    File copy_loss_del_data_tsv = CalcGeneCopyDosageSensitivity.del_data_tsv
    File copy_gain_dup_data_tsv = CalcGeneCopyDosageSensitivity.dup_data_tsv
    File? diagnostics = MergeDiagnostics.tarball
  }
}


# Filter GTF to CDSs and gene bodies
task FilterGtf {
  input {
    File gtf
    File gtf_idx

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
    zcat ~{gtf} | awk '{ if ($3 == "CDS") print $0 }' | bgzip -c \
    > ~{gtf_prefix}.coding.gtf.gz
    tabix -f ~{gtf_prefix}.coding.gtf.gz

    # Subset GTF to gene bodies
    zcat ~{gtf} | awk '{ if ($3 == "gene") print $0 }' | bgzip -c \
    > ~{gtf_prefix}.genes.gtf.gz
    tabix -f ~{gtf_prefix}.genes.gtf.gz
  >>>

  output {
    File coding_gtf = "~{gtf_prefix}.coding.gtf.gz"
    File coding_gtf_idx = "~{gtf_prefix}.coding.gtf.gz.tbi"
    File genes_gtf = "~{gtf_prefix}.genes.gtf.gz"
    File genes_gtf_idx = "~{gtf_prefix}.genes.gtf.gz.tbi"
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
