#######################
#    DSMap Project    #
#######################
#
# CleanAFInfo.wdl
#
# Bin and annotate a genome or set of chromosomes
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Structs.wdl"

task BinAndAnnotateGenome {
	input {
    File contigs_fai
    File bin_annotations_list
    File pair_annotations_list
    File bin_exclusion_mask
    String prefix

    Int bin_size
    Int max_pair_distance
    Int bins_per_shard
    Int pairs_for_pca

    String athena_docker

    RuntimeAttr? runtime_attr_make_bins
  }

  Array[Array[String]] contigs = read_tsv(contigs_fai)


  # Step 1. Create all 1D bins
  call MakeBins {
    input:
      contigs_fai=contigs_fai,
      bin_size=bin_size,
      bin_exclusion_mask=bin_exclusion_mask,
      prefix=prefix,
      athena_docker=athena_docker,
      runtime_attr_override=runtime_attr_make_bins
  }

  # Scatter over contigs
  
    # Step 2. Annotate 1D bins per contig

    # Step 3. Pair 2D bins and add 2D annotations

    # Step 4. Uniformly sample pairs per chromosome

    # (Step 5 called outside of scatter; see below)

    # Step 6. Apply PCA transformation to 2D pairs

    # Step 7. Merge & sort all 2D bins per chromosome

  # Step 5. Learn PCA transformation from sampled bins

}


# Create sequential uniform bins for all chromosomes
task MakeBins {
  input {
    File contigs_fai
    Int bin_size
    File bin_exclusion_mask
    String prefix
    
    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Create .genome file from contigs_fai
    cut -f1-2 ~{contigs_fai} > contigs.genome

    # Make bins & tabix output BED
    athena make-bins \
      --exclude-all ~{bin_exclusion_mask} \
      --buffer ~{bin_size} \
      --include-chroms $( cut -f1 contigs.genome | paste -s -d, ) \
      --bgzip \
      contigs.genome \
      ~{bin_size} \
      ~{prefix}.bins.bed.gz
    tabix -f ~{prefix}.bins.bed.gz
  }

  output {
    File bins_bed = "~{prefix}.bins.bed.gz"
    File bins_bed_idx = "~{prefix}.bins.bed.gz.tbi"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: dsmap_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
