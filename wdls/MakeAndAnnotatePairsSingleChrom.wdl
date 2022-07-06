#######################
#    DSMap Project    #
#######################
#
# MakeAndAnnotatePairsSingleChrom.wdl
#
# Parallelized creation and annotation of pairs of bins on a single chromosome
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow MakeAndAnnotatePairsSingleChrom {
  input {
    File bins
    File bins_idx
    File bedtools_genome_file
    File? pair_exclusion_mask
    String contig
    Int max_pair_distance
    Int shard_size
    String prefix

    File? pair_annotations_list_localize
    File? pair_annotations_list_remote
    File? pair_annotations_list_ucsc
    String? ref_build
    String? ref_fasta
    Int? bin_size

    Boolean sample_pairs_for_pca = true
    Int? pairs_to_sample_for_pca = 0

    String athena_docker
    String athena_cloud_docker

    RuntimeAttr? runtime_attr_chrom_shard
    RuntimeAttr? runtime_attr_make_pairs
    RuntimeAttr? runtime_attr_annotate_pairs
    RuntimeAttr? runtime_attr_sample_pairs
    RuntimeAttr? runtime_attr_merge_annotated_pairs
  }

  # Shard bins for parallelized pairing & annotation
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

  # Scatter over sharded bins for pairing & annotation
  scatter ( shard in ChromShard.shards ) {
    # Make pairs per shard
    call MakePairs {
      input:
        query_bins=shard,
        all_bins=bins,
        all_bins_idx=bins_idx,
        pair_exclusion_mask=pair_exclusion_mask,
        max_pair_distance=max_pair_distance,
        prefix=prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_make_pairs
    }

    # Annotate pairs per shard
    call AnnotatePairs {
      input:
        pairs=MakePairs.pairs,
        pairs_idx=MakePairs.pairs_idx,
        bedtools_genome_file=bedtools_genome_file,
        pair_annotations_list_localize=pair_annotations_list_localize,
        pair_annotations_list_remote=pair_annotations_list_remote,
        pair_annotations_list_ucsc=pair_annotations_list_ucsc,
        ref_build=ref_build,
        ref_fasta=ref_fasta,
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_override=runtime_attr_annotate_pairs
    }
  }

  # Merge all pairs
  call Utils.MergeBEDs as MergeAnnotatedPairs {
    input:
      beds=AnnotatePairs.annotated_pairs,
      prefix="~{prefix}.annotated_pairs.~{contig}",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_annotated_pairs
  }


  # Sample N random pairs for feature PCA, if optioned
  if ( sample_pairs_for_pca ) {
    call SamplePairs {
      input:
        annotated_pairs=MergeAnnotatedPairs.merged_bed,
        annotated_pairs_idx=MergeAnnotatedPairs.merged_bed_idx,
        sample_size=pairs_to_sample_for_pca,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_sample_pairs
    }
  }


  output {

    File annotated_pairs = MergeAnnotatedPairs.merged_bed
    File annotated_pairs_idx = MergeAnnotatedPairs.merged_bed_idx

    File? downsampled_pairs = SamplePairs.sampled_pairs
    File? downsampled_pairs_idx = SamplePairs.sampled_pairs_idx
    
  }
}


# Create pairs for an input BED file of bins
task MakePairs {
  input {
    File query_bins
    File all_bins
    File all_bins_idx
    Int max_pair_distance
    String prefix
    File? pair_exclusion_mask

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail
    touch empty.txt

    # Build options for athena pair-bins
    athena_options=""
    if [ "~{defined(pair_exclusion_mask)}" == "true" ]; then
      athena_options="$athena_options --exclusion-list ~{default='empty.txt' pair_exclusion_mask}"
    fi

    # Pair bins with athena
    athena_cmd="athena pair-bins --bgzip --max-dist ~{max_pair_distance}"
    athena_cmd="$athena_cmd --annotate-distance --sort-features --annotate-absdiff"
    athena_cmd="$athena_cmd $athena_options --bin-superset ~{all_bins}"
    athena_cmd="$athena_cmd ~{query_bins} ~{prefix}.pairs.bed.gz"
    echo -e "Now pairing bins using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ~{prefix}.pairs.bed.gz
  }

  output {
    File pairs = "~{prefix}.pairs.bed.gz"
    File pairs_idx = "~{prefix}.pairs.bed.gz.tbi"
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


# Add annotations to bin-pairs generated by MakePairs
task AnnotatePairs {
  input {
    File pairs
    File pairs_idx
    File bedtools_genome_file

    File? pair_annotations_list_localize
    File? pair_annotations_list_remote
    File? pair_annotations_list_ucsc
    String? ref_build
    File? ref_fasta
    Int? bin_size

    Int? query_slop = 1000

    String athena_cloud_docker

    RuntimeAttr? runtime_attr_override
  }

  String out_prefix = basename(pairs, ".bed.gz")
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 10 + (10 * ceil(size(pairs, "GB"))),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail
    touch empty.txt

    # Prepare inputs
    zcat ~{pairs} \
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
    if [ "~{defined(pair_annotations_list_localize)}" == "true" ]; then
      cut -f1 ~{pair_annotations_list_localize} \
      | gsutil -m cp -I ./
      while read path action tname; do
        echo -e "$( basename $path )\t$action\t$tname" >> local_tracks.pairs.tsv
      done < ~{default="empty.txt" pair_annotations_list_localize}
    fi

    # Slice large input files hosted remotely with athena slice-remote
    if [ "~{defined(pair_annotations_list_remote)}" == "true" ]; then
      if [ "~{defined(ref_fasta)}" == "true" ]; then
        remote_options="--ref-fasta ~{default='empty.txt' ref_fasta}"
      else
        remote_options=""
      fi
      athena slice-remote $remote_options \
        --updated-tsv local_slices.pairs.tsv \
        ~{default="empty.txt" pair_annotations_list_remote} \
        regions.bed
      cat local_slices.pairs.tsv >> local_tracks.pairs.tsv
    fi

    # Build options for athena annotate-pairs
    athena_options=""
    if [ -s local_tracks.pairs.tsv ]; then
      athena_options="$athena_options --track-list local_tracks.pairs.tsv"
    fi
    if [ "~{defined(pair_annotations_list_ucsc)}" == "true" ]; then
      athena_options="$athena_options --ucsc-list ~{default='empty.txt' pair_annotations_list_ucsc}"
    fi
    if [ "~{defined(ref_build)}" == "true" ]; then
      athena_options="$athena_options --ucsc-ref ~{default='' ref_build}"
    fi
    if [ "~{defined(ref_fasta)}" == "true" ]; then
      athena_options="$athena_options --fasta ~{default='empty.txt' ref_fasta}"
    fi
    if [ "~{defined(bin_size)}" == "true" ]; then
      athena_options="$athena_options --binsize ~{default='' bin_size}"
    fi

    # Annotate pairs with athena
    athena_cmd="athena annotate-pairs --bgzip --no-ucsc-chromsplit $athena_options"
    athena_cmd="$athena_cmd ~{pairs} ~{out_prefix}.annotated.pairs.bed.gz"
    echo -e "Now pairing bins using command:\n$athena_cmd"
    eval $athena_cmd 2> >(fgrep -v "is marked as paired, but its mate does not occur" >&2)
    tabix -f ~{out_prefix}.annotated.pairs.bed.gz
  }

  output {
    File annotated_pairs = "~{out_prefix}.annotated.pairs.bed.gz"
    File annotated_pairs_idx = "~{out_prefix}.annotated.pairs.bed.gz.tbi"
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


# Randomly sample a subset of pairs from a BED for PCA in parent workflow
task SamplePairs {
  input {
    File annotated_pairs
    File annotated_pairs_idx
    Int? sample_size
    
    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  String out_prefix = basename(annotated_pairs, ".bed.gz")
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 10 + (10 * ceil(size(annotated_pairs, "GB"))),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    tabix -H ~{annotated_pairs} > header.bed

    zcat ~{annotated_pairs} | grep -ve '^#' \
    | shuf -n ~{sample_size} --random-source=<( zcat ~{annotated_pairs} ) \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.bed - \
    | bgzip -c \
    > ~{out_prefix}.downsampled.bed.gz
    tabix -f ~{out_prefix}.downsampled.bed.gz
  }

  output {
    File sampled_pairs = "~{out_prefix}.downsampled.bed.gz"
    File sampled_pairs_idx = "~{out_prefix}.downsampled.bed.gz.tbi"
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
