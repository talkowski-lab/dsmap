#######################
#    DSMap Project    #
#######################
#
# CountCnvsInBins.wdl
#
# Count CNVs per 1D bin or 2D bin-pair
#
# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow CountCnvsInBins {
  input {
    # General inputs
    File contigs_fai
    String prefix

    # CNV inputs
    String del_vcf
    File del_vcf_idx
    String dup_vcf
    File dup_vcf_idx

    # Bin or bin-pair inputs
    String bins_bucket
    String bins_bed_prefix
    Boolean bins_are_paired

    # Count options
    Boolean count_probs
    Boolean full_segment_overlap = false

    # Diagnostics options
    Boolean run_diagnostics = true

    # Dockers
    String athena_docker
    String athena_cloud_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_count_bin_cnvs
    RuntimeAttr? runtime_attr_diagnostics
  }


  String bin_pair_prefix = if bins_are_paired then "pairs" else "bins"
  String count_probs_prefix = if count_probs then "probs" else "counts"
  String del_prefix = prefix + ".DEL" + "." + bin_pair_prefix + "." + count_probs_prefix
  String dup_prefix = prefix + ".DUP" + "." + bin_pair_prefix + "." + count_probs_prefix

  # If inputs are for bin-pairs, run Steps 1-3
  if ( bins_are_paired ) {

    Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

    # Parallelize over chromosomes as bin-pair input BEDs are split by chromosome
    scatter ( contig in contigs ) {

      # Step 1a. Subset DEL VCF to chromosome
      call SubsetVCFSingleChrom as SubsetDelSingleChrom {
        input:
          vcf=del_vcf,
          contig=contig,
          athena_cloud_docker=athena_cloud_docker,
          runtime_attr_override=runtime_attr_subset_vcf
      }

      # Step 1b. Subset DUP VCF to chromosome
      call SubsetVCFSingleChrom as SubsetDupSingleChrom {
        input:
          vcf=dup_vcf,
          contig=contig,
          athena_cloud_docker=athena_cloud_docker,
          runtime_attr_override=runtime_attr_subset_vcf
      }

      # Infer pairs BED file path
      File bins_bed = bins_bucket + "/" + bins_bed_prefix + "." + contig + ".bed.gz" 
      File bins_bed_idx = bins_bed + ".tbi"

      # Step 2a. Count DELs with breakpoints in 2D bin-pairs
      call CountCnvs as CountPairDels {
        input:
          vcf=SubsetDelSingleChrom.chr_vcf,
          vcf_idx=SubsetDelSingleChrom.chr_vcf_idx,
          bins_bed=bins_bed,
          bins_bed_idx=bins_bed_idx,
          bins_are_paired=bins_are_paired,
          count_probs=count_probs,
          contig=contig,
          prefix=del_prefix,
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_count_bin_cnvs
      }

      # Step 2b. Count DUPs with breakpoints in 2D bin-pairs
      call CountCnvs as CountPairDups {
        input:
          vcf=SubsetDupSingleChrom.chr_vcf,
          vcf_idx=SubsetDupSingleChrom.chr_vcf_idx,
          bins_bed=bins_bed,
          bins_bed_idx=bins_bed_idx,
          bins_are_paired=bins_are_paired,
          count_probs=count_probs,
          contig=contig,
          prefix=dup_prefix,
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_count_bin_cnvs
      }
    }

    # [Optional] Step 3. Run diagnostics
    if ( run_diagnostics ) {

      # Step 3a. Run diagnostics on DELs with breakpoints in 2D bins
      call Utils.GetPairDiagnostics as GetDelPairDiagnostics {
        input:
          pair_counts=CountPairDels.bins_w_counts,
          counts_are_probs=count_probs,
          cnv="DEL",
          prefix=del_prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 3b. Run diagnostics on DUPs with breakpoints in 2D bins
      call Utils.GetPairDiagnostics as GetDupPairDiagnostics {
        input:
          pair_counts=CountPairDups.bins_w_counts,
          counts_are_probs=count_probs,
          cnv="DUP",
          prefix=dup_prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 3c. Tar diagnostics for convenience
      call Utils.MakeTarball as MergeDelPairDiagnostics {
        input:
          files_to_tar=GetDelPairDiagnostics.outputs,
          tarball_prefix="~{prefix}.DEL.CountCnvsInBins.~{bin_pair_prefix}.~{count_probs_prefix}.diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
      call Utils.MakeTarball as MergeDupPairDiagnostics {
        input:
          files_to_tar=GetDupPairDiagnostics.outputs,
          tarball_prefix="~{prefix}.DUP.CountCnvsInBins.~{bin_pair_prefix}.~{count_probs_prefix}.diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
    }
  }

  # If inputs are for bin-pairs, run Steps 4-5
  if ( !bins_are_paired ) {

    # Infer bins BED file path
    File bins_bed = bins_bucket + "/" + bins_bed_prefix + ".bed.gz" 
    File bins_bed_idx = bins_bed + ".tbi"

    # Step 4a. Count DELs overlapping 1D bins
    call CountCnvs as CountBinDels {
      input:
        vcf=del_vcf,
        vcf_idx=del_vcf_idx,
        bins_bed=bins_bed,
        bins_bed_idx=bins_bed_idx,
        bins_are_paired=bins_are_paired,
        count_probs=count_probs,
        full_segment_overlap=full_segment_overlap,
        prefix=del_prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_bin_cnvs
    }

    # Step 4b. Count DUPs overlapping 1D bins
    call CountCnvs as CountBinDups {
      input:
        vcf=dup_vcf,
        vcf_idx=dup_vcf_idx,
        bins_bed=bins_bed,
        bins_bed_idx=bins_bed_idx,
        bins_are_paired=bins_are_paired,
        count_probs=count_probs,
        full_segment_overlap=full_segment_overlap,
        prefix=dup_prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_bin_cnvs
    }

    # [Optional] Step 5. Run diagnostics
    if ( run_diagnostics ) {

      # Step 5a. Run diagnostics on DELs overlapping 1D bins
      call Utils.GetBinDiagnostics as GetDelBinDiagnostics {
        input:
          bin_counts=CountBinDels.bins_w_counts,
          cnv="DEL",
          prefix=del_prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 5b. Run diagnostics on DUPs overlapping 1D bins
      call Utils.GetBinDiagnostics as GetDupBinDiagnostics {
        input:
          bin_counts=CountBinDups.bins_w_counts,
          cnv="DUP",
          prefix=dup_prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 5c. Tar diagnostics for convenience
      call Utils.MakeTarball as MergeDelBinDiagnostics {
        input:
          files_to_tar=GetDelBinDiagnostics.outputs,
          tarball_prefix="~{prefix}.DEL.CountCnvsInBins.~{bin_pair_prefix}.~{count_probs_prefix}.diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
      call Utils.MakeTarball as MergeDupBinDiagnostics {
        input:
          files_to_tar=GetDupBinDiagnostics.outputs,
          tarball_prefix="~{prefix}.DUP.CountCnvsInBins.~{bin_pair_prefix}.~{count_probs_prefix}.diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
    }
  }


  output {
    # Resolve optional types
    Array[File] bin_del_counts = select_first(
      [CountPairDels.bins_w_counts, [select_first([CountBinDels.bins_w_counts, ""])]]
    )
    Array[File] bin_del_counts_idxs = select_first(
      [CountPairDels.bins_w_counts_idx, [select_first([CountBinDels.bins_w_counts_idx, ""])]]
    )
    Array[File] bin_dup_counts = select_first(
      [CountPairDups.bins_w_counts, [select_first([CountBinDups.bins_w_counts, ""])]]
    )
    Array[File] bin_dup_counts_idxs = select_first(
      [CountPairDups.bins_w_counts_idx, [select_first([CountBinDups.bins_w_counts_idx, ""])]]
    )
    
    File? bin_del_diagnostics = (
      if bins_are_paired then MergeDelPairDiagnostics.tarball
      else MergeDelBinDiagnostics.tarball
    )
    File? bin_dup_diagnostics = (
      if bins_are_paired then MergeDupPairDiagnostics.tarball
      else MergeDupBinDiagnostics.tarball
    )
  }
}


# Subset VCF to a single chromosome
task SubsetVCFSingleChrom {
  input {
    String vcf  # VCF passed as string to allow for remote tabixing without localizing entire VCF
    String contig

    String athena_cloud_docker

    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf, ".vcf.gz")

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size(vcf, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Localize variants from contig
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{prefix}.~{contig}.vcf.gz
    tabix -f ~{prefix}.~{contig}.vcf.gz

  >>>

  output {
    File chr_vcf = "~{prefix}.~{contig}.vcf.gz"
    File chr_vcf_idx = "~{prefix}.~{contig}.vcf.gz.tbi"
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


# Count number of qualifying CNVs per segment
task CountCnvs {
  input {
    File vcf
    File vcf_idx
    File bins_bed
    File bins_bed_idx
    Boolean bins_are_paired
    Boolean count_probs
    Boolean full_segment_overlap = false
    String? contig
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  String contig_prefix = select_first([contig, ""])
  String outfile = (
    prefix +
    ( if contig_prefix == "" then contig_prefix else "." + contig_prefix ) +
    ".bed.gz"
  )

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10 + ceil(2 * size([bins_bed, vcf], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Create option for filtering to CNVs with full bin overlap
    athena_options=""
    if [ "~{full_segment_overlap}" == "true" ] && [ "~{bins_are_paired}" == "false" ]; then
      athena_options="$athena_options --fraction 1.0"
    fi

    # Count SVs
    athena_cmd="athena count-sv --query-format ~{true='pairs' false='bins' bins_are_paired}"
    athena_cmd="$athena_cmd --comparison ~{true='breakpoint' false='overlap' bins_are_paired}"
    athena_cmd="$athena_cmd ~{true='--probabilities' false='' count_probs}"
    athena_cmd="$athena_cmd $athena_options"
    athena_cmd="$athena_cmd --outfile ~{outfile} --bgzip"
    athena_cmd="$athena_cmd ~{vcf} ~{bins_bed}"
    echo -e "Now counting SVs using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ~{outfile}
  >>>

  output {
    File bins_w_counts = outfile
    File bins_w_counts_idx = outfile + ".tbi"
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
