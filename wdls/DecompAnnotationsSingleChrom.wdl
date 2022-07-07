#######################
#    DSMap Project    #
#######################
#
# DecompAnnotationsSingleChrom.wdl
#
# Parallelized decomposition of annotations for all bin-pairs on a single chromosome
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow DecompAnnotationsSingleChrom {
  input {
    File pairs
    File pairs_idx
    File pca_model
    String contig
    Int shard_size = 25000
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_chrom_shard
    RuntimeAttr? runtime_attr_apply_pca
    RuntimeAttr? runtime_attr_merge_pairs
  }

  # Shard pairs for parallelized decomposition
  call Utils.SingleChromShard as ChromShard {
    input:
      infile=pairs,
      infile_idx=pairs_idx,
      contig=contig,
      shard_size=shard_size,
      prefix=prefix,
      file_format="bed",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_chrom_shard
  }

  # Apply pre-learned PCA model to each shard
  scatter ( shard in ChromShard.shards ) {
    call DecompAnnos {
      input:
        bed=shard,
        pca_model=pca_model,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_apply_pca
    }
  }

  # Merge decomposed shards
  call Utils.MergeBEDs as MergeDecompedPairs {
    input:
      beds=DecompAnnos.decomped_bed,
      prefix="~{prefix}.eigen.~{contig}",
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_merge_pairs
  }

  output {
    File decomped_bed = MergeDecompedPairs.merged_bed
    File decomped_bed_idx = MergeDecompedPairs.merged_bed_idx
  }
}


# Apply a pre-trained PCA model to an input BED
task DecompAnnos {
  input {
    File bed
    File pca_model

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  String out_prefix = basename(bed, ".bed.gz")
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 10 + (10 * ceil(size(bed, "GB"))),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -eu -o pipefail

    # Apply PCA model to annotations with athena
    athena_cmd="athena eigen-bins --bgzip"
    athena_cmd="$athena_cmd -o ~{out_prefix}.eigen.bed.gz"
    athena_cmd="$athena_cmd --precomputed-parameters ~{pca_model}"
    athena_cmd="$athena_cmd ~{bed}"
    echo -e "Now applying PCA decomposition using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f "~{out_prefix}.eigen.bed.gz"
  }

  output {
    File decomped_bed = "~{out_prefix}.eigen.bed.gz"
    File decomped_bed_idx = "~{out_prefix}.eigen.bed.gz.tbi"
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
