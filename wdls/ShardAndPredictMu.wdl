#######################
#    DSMap Project    #
#######################
#
# ShardAndPredictMu.wdl
#
# Predict SV mutation rates in parallel by sharding input
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow ShardAndPredictMu {
  input {
    File bed
    File bed_idx
    File trained_model
    Int shard_size
    String contig
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_shard_bed
    RuntimeAttr? runtime_attr_apply_mu_model
    RuntimeAttr? runtime_attr_merge_beds
  }

  # Shard input bed by shard_size
  call Utils.SingleChromShard as ShardBED {
    input:
      infile=bed,
      infile_idx=bed_idx,
      contig=contig,
      shard_size=shard_size,
      prefix=prefix,
      file_format="bed",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_shard_bed
  }

  # Apply model per shard
  scatter ( shard in ShardBED.shards ) {
    call ApplyMuModel {
      input:
        pairs_bed=shard,
        trained_model=trained_model,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_apply_mu_model
    }
  }

  # Merge shards
  call Utils.MergeBEDs as MergeBEDs {
    input:
      beds=ApplyMuModel.bed_w_mu,
      prefix="~{prefix}.~{contig}.mu",
      beds_are_bgzipped=true,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_beds
  }

  output {
    File pairs_w_mu = MergeBEDs.merged_bed
    File pairs_w_mu_idx = MergeBEDs.merged_bed_idx
  }
}


# Apply a pre-trained mutation rate model on a set of input bin-pairs
task ApplyMuModel {
  input {
    File pairs_bed
    File trained_model

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  String out_prefix = basename(pairs_bed, '.bed.gz')

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size(pairs_bed, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Predict mutation rates for all bin-pairs
    athena mu-predict \
      --trained-model ~{trained_model} \
      --outfile ~{out_prefix}.mu.bed.gz \
      --bgzip \
      ~{pairs_bed}
    tabix -f ~{out_prefix}.mu.bed.gz
  }

  output {
    File bed_w_mu = "~{out_prefix}.mu.bed.gz"
    File bed_w_mu_idx = "~{out_prefix}.mu.bed.gz.tbi"
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

