#######################
#    DSMap Project    #
#######################
#
# AnnotateBinsSingleChrom.wdl
#
# Parallelized annotation of  bins on a single chromosome
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow AnnotateBinsSingleChrom {
  input {
    File bins
    File bins_idx
    File bedtools_genome_file
    String contig
    Int shard_size
    String prefix

    File? bin_annotations_list_localize
    File? bin_annotations_list_remote
    File? bin_annotations_list_ucsc
    String? ref_build = "hg38"
    File? ref_fasta
    File? snv_mutrates_tsv

    String athena_docker
    String athena_cloud_docker

    RuntimeAttr? runtime_attr_chrom_shard
    RuntimeAttr? runtime_attr_annotate_bins
    RuntimeAttr? runtime_attr_merge_annotated_bins
  }

  # Shard bins for parallelized annotation
  call Utils.SingleChromShard as ChromShard {
    input:
      infile=bins,
      infile_idx=bins_idx,
      contig=contig,
      shard_size=shard_size,
      prefix=prefix,
      file_format="bed",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_chrom_shard
  }

  # Annotate each shard
  scatter ( shard in ChromShard.shards ) {
    call AnnotateBins {
      input:
        bed=shard,
        bedtools_genome_file=bedtools_genome_file,
        bin_annotations_list_localize=bin_annotations_list_localize,
        bin_annotations_list_remote=bin_annotations_list_remote,
        bin_annotations_list_ucsc=bin_annotations_list_ucsc,
        ref_build=ref_build,
        ref_fasta=ref_fasta,
        snv_mutrates_tsv=snv_mutrates_tsv,
        athena_cloud_docker=athena_cloud_docker
    }
  }

  # Merge annotated shards
  call Utils.MergeBEDs as MergeAnnotatedBins {
    input:
      beds=AnnotateBins.annotated_bed,
      prefix="~{prefix}.annotated_bins.~{contig}",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_annotated_bins
  }

  output {
    File annotated_bins = MergeAnnotatedBins.merged_bed
    File annotated_bins_idx = MergeAnnotatedBins.merged_bed_idx
  }
}


# Annotate a BED file of  bins
task AnnotateBins {
  input {
    File bed
    File bedtools_genome_file
    
    File? bin_annotations_list_localize
    File? bin_annotations_list_remote
    File? bin_annotations_list_ucsc
    File? ref_fasta
    String? ref_build
    File? snv_mutrates_tsv

    Int? query_slop = 1000

    String athena_cloud_docker

    RuntimeAttr? runtime_attr_override
  }

  String out_prefix = basename(bed, ".bed.gz")
  
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
    touch empty.txt

    # Prepare inputs
    zcat ~{bed} \
    | fgrep -v "#" \
    | bedtools slop -i - -g ~{bedtools_genome_file} -b ~{query_slop} \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    > regions.bed
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    if [ "~{defined(ref_fasta)}" == "true" ]; then
      samtools faidx ~{default="empty.txt" ref_fasta}
    fi

    # Localize small input files hosted remotely
    if [ "~{defined(bin_annotations_list_localize)}" == "true" ]; then
      cut -f1 ~{bin_annotations_list_localize} \
      | gsutil -m cp -I ./
      while read path action tname; do
        echo -e "$( basename $path )\t$action\t$tname" >> local_tracks.tsv
      done < ~{default="empty.txt" bin_annotations_list_localize}
    fi

    # Slice large input files hosted remotely with athena slice-remote
    if [ "~{defined(bin_annotations_list_remote)}" == "true" ]; then
      if [ "~{defined(ref_fasta)}" == "true" ]; then
        remote_options="--ref-fasta ~{default='empty.txt' ref_fasta}"
      else
        remote_options=""
      fi
      athena slice-remote $remote_options \
        --updated-tsv local_slices.tsv \
        ~{default="empty.txt" bin_annotations_list_remote} \
        regions.bed
      cat local_slices.tsv >> local_tracks.tsv
    fi

    # Build options for athena annotate-bins
    athena_options=""
    if [ ! -z local_tracks.tsv ]; then
      athena_options="$athena_options --track-list local_tracks.tsv"
    fi
    if [ "~{defined(bin_annotations_list_ucsc)}" == "true" ]; then
      athena_options="$athena_options --ucsc-list ~{default='empty.txt' bin_annotations_list_ucsc}"
    fi
    if [ "~{defined(ref_fasta)}" == "true" ]; then
      athena_options="$athena_options --fasta ~{default='empty.txt' ref_fasta}"
    fi
    if [ "~{defined(snv_mutrates_tsv)}" == "true" ]; then
      athena_options="$athena_options --snv-mutrate ~{default='empty.txt' snv_mutrates_tsv}"
    fi

    # Annotate bins with athena
    athena_cmd="athena annotate-bins --ucsc-ref ~{ref_build} $athena_options"
    athena_cmd="$athena_cmd --no-ucsc-chromsplit --bgzip"
    athena_cmd="$athena_cmd ~{bed} ~{out_prefix}.annotated.bed.gz"
    echo -e "Now annotating using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ~{out_prefix}.annotated.bed.gz
  }

  output {
    File annotated_bed = "~{out_prefix}.annotated.bed.gz"
    File annotated_bed_idx = "~{out_prefix}.annotated.bed.gz.tbi"
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
