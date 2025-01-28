#######################
#    DSMap Project    #
#######################
#
# CalcSegmentDosageSensitivity.wdl
#
# Compute dosage sensitivity statistics for genomic segments in a single SV dataset
#
# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>


version 1.0

import "CountCnvsInBins.wdl"
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
    String athena_cloud_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_filter_query_chrom
    RuntimeAttr? runtime_attr_query_mu
    RuntimeAttr? runtime_attr_expand_query
    RuntimeAttr? runtime_attr_count_cnvs
    RuntimeAttr? runtime_attr_merge_data
    RuntimeAttr? runtime_attr_plot_mu_hist
    RuntimeAttr? runtime_attr_merge_diagnostics
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Create athena options
  if ( full_segment_overlap ) {
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

    # Step 2a. Tally mutation rates for all deletions in bin-pairs overlapping
    # bins housing each segment
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

    # Step 2b. Tally mutation rates for all duplications in bin-pairs overlapping
    # bins housing each segment
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

    # Step 3. To match mu query and CNV counting strategy, expand query segment
    # boundaries to nearest neighboring larger bins in mu matrix
    # NOTE: Segments with the same such neighboring bin boundaries will have
    # same mu, CNV count, and O/E estimates
    # NOTE: Assumes that mu matrices for deletions and duplications are defined
    # over the same space
    call ExpandQueryToBins {
      input:
        query=FilterQuerySingleChrom.query_chrom,
        query_idx=FilterQuerySingleChrom.query_chrom_idx,
        mu_bed=del_mu_bed,
        mu_bed_idx=del_mu_bed_idx,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_expand_query
    }

    # Step 4. Count deletions and duplications overlapping bins housing each segment
    call CountCnvsInBins.CountCnvsInBins as CountQueryCnvs {
      input:
        del_vcf=del_vcf,
        del_vcf_idx=del_vcf_idx,
        dup_vcf=dup_vcf,
        dup_vcf_idx=dup_vcf_idx,
        bins_bucket=ExpandQueryToBins.expanded_query,
        bins_bed_prefix=basename(ExpandQueryToBins.expanded_query, ".bed.gz"),
        bins_are_paired=false,
        contigs_fai=contigs_fai,
        prefix=prefix + "." + contig,
        count_probs=false,
        full_segment_overlap=full_segment_overlap,
        run_diagnostics=false,
        athena_docker=athena_docker,
        athena_cloud_docker=athena_cloud_docker,
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_count_bin_cnvs=runtime_attr_count_cnvs
    }
  }

  # Step 5a. Merge and analyze deletion outputs
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeDelData {
    input:
      mu_tsvs=QueryMuDel.mu_tsv,
      counts_tsvs=flatten(CountQueryCnvs.bin_del_counts),
      prefix=prefix + del_prefix,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }
  
  # Step 5b. Merge and analyze duplication outputs
  # Note: for now, just merge & joint outputs. TODO: add analysis components
  call MergeMuAndCounts as MergeDupData {
    input:
      mu_tsvs=QueryMuDup.mu_tsv,
      counts_tsvs=flatten(CountQueryCnvs.bin_dup_counts),
      prefix=prefix + dup_prefix,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_data
  }

  # Gather diagnostics, if optioned
  if ( run_diagnostics ) {

    # Plot mutation rates per segment
    call Utils.PlotMuHist as PlotDelMu {
      input:
        mu_tsv=MergeDelData.merged_tsv,
        cnv="DEL",
        x_title="Deletions per allele per generation",
        y_title="Segments",
        out_prefix=prefix + del_prefix,
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_plot_mu_hist
    }
    call Utils.PlotMuHist as PlotDupMu {
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
        files_to_tar=[PlotDelMu.mu_hist, PlotDupMu.mu_hist],
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

    # Sort by chr, start, end position
    zcat ~{query} | grep -e "^#" | sed -n '1p' > query.header
    tabix -h ~{query} ~{contig} | grep -ve "^#" | sort -Vk1,1 -k2,2n -k3,3n \
    | cat query.header - \
    | bgzip -c > ~{query_prefix}.~{contig}.bed.gz
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


# To match mu and CNV counting strategies, extend query segments to boundaries of
# nearest neighboring larger bins in mu matrix
# TODO: Filter out SVs from VCF that lie in bin pairs not in the mu matrix
# e.g. SVs >500kb when mu matrix is limited to <=500kb
# TODO: What if query segment is larger than the largest allowed bin pair?
task ExpandQueryToBins {
  input {
    File query
    File query_idx
    File mu_bed
    File mu_bed_idx

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  # TODO: Adjust this for other file types
  String query_prefix = basename(query, ".bed.gz")

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: 10 + ceil(2 * size([query, mu_bed], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Collect query left and right boundaries as single bp intervals
    zcat ~{query} | grep -ve "^#" | awk -v OFS="\t" '{ print $1, $2, $2+1 }' \
    | bgzip -c > query.left.bed.gz
    zcat ~{query} | grep -ve "^#" | awk -v OFS="\t" '{ print $1, $3, $3+1 }' \
    | bgzip -c > query.right.bed.gz

    # Collect all mutation rate matrix bin pair left and right boundaries
    # together as single bp intervals
    {
      ( zcat ~{mu_bed} | grep -ve "^#" | awk -v OFS="\t" '{ print $1, $2, $2+1 }' );
      ( zcat ~{mu_bed} | grep -ve "^#" | awk -v OFS="\t" '{ print $1, $3, $3+1 }' );
    } \
    | sort -Vk1,1 -k2,2n | uniq | bgzip -c > mu.bin_bounds.bed.gz

    # Expand query segments to bin pair boundaries to match mutation rate querying
    paste \
    <(
      bedtools closest -a query.left.bed.gz -b mu.bin_bounds.bed.gz -id -D ref \
      | cut -f4,5
    ) <(
      bedtools closest -a query.right.bed.gz -b mu.bin_bounds.bed.gz -iu -D ref \
      | cut -f5
    ) \
    | cat <( zcat ~{query} | sed -n '1p' ) - | bgzip -c > ~{query_prefix}.expanded.bed.gz
    tabix -f ~{query_prefix}.expanded.bed.gz
  >>>

  output {
    File expanded_query = "~{query_prefix}.expanded.bed.gz"
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
