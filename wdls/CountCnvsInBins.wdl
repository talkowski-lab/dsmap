#######################
#    DSMap Project    #
#######################
#
# CountCnvsInBins.wdl
#
# Count CNVs per 1D bin and 2D bin-pair
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
    String? bins_bed_prefix
    Boolean bins_are_paired

    # Count options
    Boolean count_probs

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


  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
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

    # Infer bins or pairs BED file path
    File bins_bed = bins_bucket + "/" + bins_bed_prefix + "." + contig + ".bed.gz" 
    File bins_bed_idx = bins_bed + ".tbi"

    # Step 2a. Count DELs overlapping 1D bins or with breakpoints in 2D bin-pairs
    call CountCnvs as CountBinDels {
      input:
        vcf=SubsetDelSingleChrom.vcf,
        vcf_idx=SubsetDelSingleChrom.vcf_idx,
        bins_bed=bins_bed,
        bins_bed_idx=bins_bed_idx,
        bins_are_paired=bins_are_paired,
        count_probs=count_probs,
        contig=contig,
        prefix="~{prefix}.DEL",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_bin_cnvs
    }

    # Step 2b. Count DUPs overlapping 1D bins or with breakpoints in 2D bin-pairs
    call CountCnvs as CountBinDups {
      input:
        vcf=SubsetDupSingleChrom.vcf,
        vcf_idx=SubsetDupSingleChrom.vcf_idx,
        bins_bed=bins_bed,
        bins_bed_idx=bins_bed_idx,
        bins_are_paired=bins_are_paired,
        count_probs=count_probs,
        contig=contig,
        prefix="~{prefix}.DUP",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_count_bin_cnvs
    }
  }

  # [Optional] Step 3. Run diagnostics
  # If 1D bins were passed in, run steps 3ai-3aiii
  # Otherwise if 2D bin-pairs were passed in, run steps 3bi-3biii
  if ( run_diagnostics ) {

    if ( !bins_are_paired ) {

      # Step 3ai. Run diagnostics on DELs overlapping 1D bins
      call Utils.GetBinDiagnostics as GetDelBinDiagnostics {
        input:
          bin_counts=CountBinDels.bins_w_counts,
          counts_are_probs=count_probs,
          cnv="DEL",
          prefix=prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 3aii. Run diagnostics on DUPs overlapping 1D bins
      call Utils.GetBinDiagnostics as GetDupBinDiagnostics {
        input:
          bin_counts=CountBinDups.bins_w_counts,
          counts_are_probs=count_probs,
          cnv="DUP",
          prefix=prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 3aiii. Tar diagnostics for convenience
      call Utils.MakeTarball as MergeDelBinDiagnostics {
        input:
          files_to_tar=GetDelBinDiagnostics.outputs,
          tarball_prefix="~{prefix}.DEL.CountCnvsInBins.bin_diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
      call Utils.MakeTarball as MergeDupBinDiagnostics {
        input:
          files_to_tar=GetDupBinDiagnostics.outputs,
          tarball_prefix="~{prefix}.DUP.CountCnvsInBins.bin_diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
    }

    if ( bins_are_paired ) {

      # Step 3bi. Run diagnostics on DELs with breakpoints in 2D bins
      call Utils.GetPairDiagnostics as GetDelPairDiagnostics {
        input:
          pair_counts=CountBinDels.bins_w_counts,
          counts_are_probs=count_probs,
          cnv="DEL",
          prefix=prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 3bii. Run diagnostics on DUPs with breakpoints in 2D bins
      call Utils.GetPairDiagnostics as GetDupPairDiagnostics {
        input:
          pair_counts=CountBinDups.bins_w_counts,
          counts_are_probs=count_probs,
          cnv="DUP",
          prefix=prefix,
          dsmap_r_docker=dsmap_r_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }

      # Step 3biii. Tar diagnostics for convenience
      call Utils.MakeTarball as MergeDelPairDiagnostics {
        input:
          files_to_tar=GetDelPairDiagnostics.outputs,
          tarball_prefix="~{prefix}.DEL.CountCnvsInBins.pair_diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
      call Utils.MakeTarball as MergeDupPairDiagnostics {
        input:
          files_to_tar=GetDupPairDiagnostics.outputs,
          tarball_prefix="~{prefix}.DUP.CountCnvsInBins.pair_diagnostics",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_diagnostics
      }
    }
  }

  output {
    Array[File] bin_del_counts = CountBinDels.bins_w_counts
    Array[File] bin_del_counts_idxs = CountBinDels.bins_w_counts_idx
    Array[File] bin_dup_counts = CountBinDups.bins_w_counts
    Array[File] bin_dup_counts_idxs = CountBinDups.bins_w_counts_idx
    File? bin_del_diagnostics = MergeDelBinDiagnostics.tarball
    File? bin_dup_diagnostics = MergeDupBinDiagnostics.tarball
    File? pair_del_diagnostics = MergeDelPairDiagnostics.tarball
    File? pair_dup_diagnostics = MergeDupPairDiagnostics.tarball
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
    File vcf = "~{prefix}.~{contig}.vcf.gz"
    File vcf_idx = "~{prefix}.~{contig}.vcf.gz.tbi"
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
    File bins_bed
    File bins_bed_idx
    Boolean bins_are_paired
    Boolean count_probs
    String contig
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  String outfile = (
    prefix + ".~{true='pairs' false='bins' bins_are_paired}" +
    ".~{true='probs' false='counts' count_probs}" + "." + contig + ".bed.gz"
  )

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size([bins_bed, vcf], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Count SVs
    athena_cmd="athena count-sv --query-format ~{true='pairs' false='bins' bins_are_paired}"
    athena_cmd="$athena_cmd --comparison ~{true='breakpoint' false='overlap' bins_are_paired}"
    athena_cmd="$athena_cmd ~{true='--probabilities' false='' count_probs}"
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
