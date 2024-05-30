#######################
#    DSMap Project    #
#######################
#
# CalcSegmentDosageSensitivity.wdl
#
# Compute dosage sensitivity statistics for genomic segments in a single SV dataset
#
# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow CalcSegmentDosageSensitivity {
  input {
    # General inputs
    File del_vcf
    File del_vcf_idx
    File dup_vcf
    File dup_vcf_idx
    File query
    File query_idx
    String mu_bucket
    String mu_bed_prefix
    File contigs_fai
    String prefix
    Boolean full_segment_overlap = false

    # Diagnostics options
    Boolean run_diagnostics = true

    # Dockers
    String athena_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_filter_query_chrom
    RuntimeAttr? runtime_attr_query_mu
    RuntimeAttr? runtime_attr_count_cnvs
    RuntimeAttr? runtime_attr_merge_data
    RuntimeAttr? runtime_attr_plot_mu_hist
    RuntimeAttr? runtime_attr_merge_diagnostics
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Create athena options
  if (full_segment_overlap) {
    String athena_full_overlap_option = "--fraction 1.0"
  }
  Array[String] athena_sv_options = select_all([athena_full_overlap_option])

  # Create file prefixes
  String del_overlap_prefix = if full_segment_overlap then ".CL" else ""
  String dup_overlap_prefix = if full_segment_overlap then ".CG" else ""
  String del_prefix = ".DEL" + del_overlap_prefix
  String dup_prefix = ".DUP" + dup_overlap_prefix

  # Parallelize per chromosome
  scatter ( contig in contigs ) {

    # Step 1. Filter query to chromosome
    call FilterQuerySingleChrom {
      input:
        query=query,
        query_idx=query_idx,
        contig=contig,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_filter_query_chrom
    }

    File del_mu_bed = mu_bucket + "/" + mu_bed_prefix + ".DEL." + contig + ".mu.bed.gz"
    File del_mu_bed_idx = del_mu_bed + ".tbi"
    File dup_mu_bed = mu_bucket + "/" + mu_bed_prefix + ".DUP." + contig + ".mu.bed.gz"
    File dup_mu_bed_idx = dup_mu_bed + ".tbi" 

    # Step 2a. Compute mutation rates for all deletions overlapping each segment
    call QueryMu as QueryMuDel {
      input:
        query=FilterQuerySingleChrom.query_chrom,
        query_idx=FilterQuerySingleChrom.query_chrom_idx,
        mu_bed=del_mu_bed,
        mu_bed_idx=del_mu_bed_idx,
        athena_query_options=athena_sv_options,
        prefix=basename(FilterQuerySingleChrom.query_chrom, ".bed.gz") + del_prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 2b. Compute mutation rates for all duplications overlapping each segment
    call QueryMu as QueryMuDup {
      input:
        query=FilterQuerySingleChrom.query_chrom,
        query_idx=FilterQuerySingleChrom.query_chrom_idx,
        mu_bed=dup_mu_bed,
        mu_bed_idx=dup_mu_bed_idx,
        athena_query_options=athena_sv_options,
        prefix=basename(FilterQuerySingleChrom.query_chrom, ".bed.gz") + dup_prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_query_mu
    }

    # Step 3a. Count deletions overlapping each segment
    call CountCnvs as CountDel {
      input:
        vcf=del_vcf,
        vcf_idx=del_vcf_idx,
        query=FilterQuerySingleChrom.query_chrom,
        query_idx=FilterQuerySingleChrom.query_chrom_idx,
        athena_countsv_options=athena_sv_options,
        prefix=basename(FilterQuerySingleChrom.query_chrom, ".bed.gz") + del_prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }

    # Step 3b. Count duplications overlapping each segment
    call CountCnvs as CountDup {
      input:
        vcf=dup_vcf,
        vcf_idx=dup_vcf_idx,
        query=FilterQuerySingleChrom.query_chrom,
        query_idx=FilterQuerySingleChrom.query_chrom_idx,
        athena_countsv_options=athena_sv_options,
        prefix=basename(FilterQuerySingleChrom.query_chrom, ".bed.gz") + dup_prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_cnvs
    }
  }

  # Step 4a. Merge and analyze outputs from 2a & 3a
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeDelData {
    input:
      mu_tsvs=QueryMuDel.mu_tsv,
      counts_tsvs=CountDel.counts_tsv,
      prefix=prefix + del_prefix,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }
  
  # Step 4b. Merge and analyze outputs from 2b & 3b
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeDupData {
    input:
      mu_tsvs=QueryMuDup.mu_tsv,
      counts_tsvs=CountDup.counts_tsv,
      prefix=prefix + dup_prefix,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }

  # Gather diagnostics, if optioned
  if ( run_diagnostics ) {

    # Plot mutation rates per segment
    call Utils.PlotMuDist as PlotDelMu {
      input:
        mu_tsv=MergeDelData.merged_tsv,
        cnv="DEL",
        x_title="Deletions per allele per generation",
        y_title="Segments",
        out_prefix=prefix + del_prefix,
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuDist as PlotDupMu {
      input:
        mu_tsv=MergeDupData.merged_tsv,
        cnv="DUP",
        x_title="Duplications per allele per generation",
        y_title="Segments",
        out_prefix=prefix + dup_prefix,
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }

    # Tar diagnostics, for convenience
    call Utils.MakeTarball as MergeDiagnostics {
      input:
        files_to_tar=[PlotDelMu.mu_dist, PlotDupMu.mu_dist],
        tarball_prefix="~{prefix}.CalcSegmentDosageSensitivity.diagnostics",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_merge_diagnostics
    }
  }

  output {
    File del_data_tsv = MergeDelData.merged_tsv
    File dup_data_tsv = MergeDupData.merged_tsv
    File? diagnostics = MergeDiagnostics.tarball
  }
}


# Filter query to a single chromosome
task FilterQuerySingleChrom {
  input {
    File query
    File query_idx
    String contig

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  # TODO: Adjust this for other file types
  String query_prefix = basename(query, ".bed.gz")

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size([query], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    tabix -h ~{query} ~{contig} | bgzip -c > ~{query_prefix}.~{contig}.bed.gz
    tabix -f ~{query_prefix}.~{contig}.bed.gz
  >>>

  output {
    File query_chrom = "~{query_prefix}.~{contig}.bed.gz"
    File query_chrom_idx = "~{query_prefix}.~{contig}.bed.gz.tbi"
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


# Compute mutation rates for an input query set of segments
task QueryMu {
  input {
    File query
    File query_idx
    File mu_bed
    File mu_bed_idx
    Array[String] athena_query_options
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 10 + ceil(2 * size([query, mu_bed], "GB")),
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
    athena_cmd="$athena_cmd ~{mu_bed} ~{query}"
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


# Count number of qualifying CNVs per segment
task CountCnvs {
  input {
    File vcf
    File vcf_idx
    File query
    File query_idx
    Array[String] athena_countsv_options
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size([query, vcf], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Count SVs
    athena_cmd="athena count-sv --query-format bed --outfile ~{prefix}.counts.tsv.gz"
    athena_cmd="$athena_cmd --bgzip  ~{sep=' ' athena_countsv_options}"
    athena_cmd="$athena_cmd ~{vcf} ~{query}"
    echo -e "Now counting SVs using command:\n$athena_cmd"
    eval $athena_cmd
  >>>

  output {
    File counts_tsv = "~{prefix}.counts.tsv.gz"
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


# Merge segment mutation rates and CNV counts
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
