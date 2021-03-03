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


workflow MakeAndAnnotatePairsSingleChrom {
  input {
    File bins
    File bins_idx
    File? pair_exclusion_mask
    String contig
    Int max_pair_distance
    Int shard_size
    String prefix
    Boolean sample_pairs_for_pca = true

    String athena_docker

    RuntimeAttr? runtime_attr_chrom_shard
    RuntimeAttr? runtime_attr_make_pairs
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
        athena_docker=athena_docker
    }

    # Annotate pairs per shard

    # Sample N random pairs per shard for feature PCA, if optioned
  }

  # Merge all pairs

  # Merge sampled pairs per shard for PCA, if optioned

}


task MakePairs {
  input {
    File query_bins
    File all_bins
    File all_bins_idx
    Int max_pair_distance
    String prefix
    File? pair_exclusion_mask

    String athena_docker

    RuntimeAttr? runtime_attr_make_pairs
  }
  
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
    docker: athena_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

