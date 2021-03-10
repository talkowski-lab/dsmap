#######################
#    DSMap Project    #
#######################
#
# BinAndAnnotateGenome.wdl
#
# Bin and annotate a genome or set of chromosomes
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"
import "AnnotateBinsSingleChrom.wdl"
import "MakeAndAnnotatePairsSingleChrom.wdl"


workflow BinAndAnnotateGenome {
  input {
    # General inputs
    File contigs_fai
    String prefix
    String? ref_build = "hg38"
    File? ref_fasta
    Boolean? visualize_features = false

    # Bin inputs
    File bin_exclusion_mask
    Int bin_size
    Int bins_per_shard
    File? bin_annotations_list_localize
    File? bin_annotations_list_remote
    File? bin_annotations_list_ucsc
    File? snv_mutrates_tsv
    
    # Pair inputs
    File? pair_exclusion_mask
    Int max_pair_distance
    Int bins_per_pair_shard
    File? pair_annotations_list_localize
    File? pair_annotations_list_remote
    File? pair_annotations_list_ucsc
    # Int pairs_for_pca

    # Dockers
    String athena_docker
    String athena_cloud_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_make_bins
    RuntimeAttr? runtime_attr_chrom_shard
    RuntimeAttr? runtime_attr_annotate_bins
    RuntimeAttr? runtime_attr_merge_annotated_bins
    RuntimeAttr? runtime_attr_make_pairs
    RuntimeAttr? runtime_attr_annotate_pairs
    RuntimeAttr? runtime_attr_merge_annotated_pairs
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

  # Process each chromosome in parallel
  scatter ( contig in contigs ) {

    # Step 2. Annotate 1D bins per contig
    call AnnotateBinsSingleChrom.AnnotateBinsSingleChrom as AnnotateBins {
      input:
        bins=MakeBins.bins_bed,
        bins_idx=MakeBins.bins_bed_idx,
        bedtools_genome_file=MakeBins.bedtools_genome_file,
        contig=contig[0],
        shard_size=bins_per_shard,
        prefix=prefix,
        bin_annotations_list_localize=bin_annotations_list_localize,
        bin_annotations_list_remote=bin_annotations_list_remote,
        bin_annotations_list_ucsc=bin_annotations_list_ucsc,
        ref_build=ref_build,
        ref_fasta=ref_fasta,
        snv_mutrates_tsv=snv_mutrates_tsv,
        athena_docker=athena_docker,
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_chrom_shard=runtime_attr_chrom_shard,
        runtime_attr_annotate_bins=runtime_attr_annotate_bins,
        runtime_attr_merge_annotated_bins=runtime_attr_merge_annotated_bins
    }

    # Step 3. Pair 2D bins and add 2D annotations
    call MakeAndAnnotatePairsSingleChrom.MakeAndAnnotatePairsSingleChrom as MakeAndAnnotatePairs {
      input:
        bins=AnnotateBins.annotated_bins,
        bins_idx=AnnotateBins.annotated_bins_idx,
        bedtools_genome_file=MakeBins.bedtools_genome_file,
        pair_exclusion_mask=pair_exclusion_mask,
        contig=contig[0],
        max_pair_distance=max_pair_distance,
        shard_size=bins_per_pair_shard,
        prefix=prefix,
        pair_annotations_list_localize=pair_annotations_list_localize,
        pair_annotations_list_remote=pair_annotations_list_remote,
        pair_annotations_list_ucsc=pair_annotations_list_ucsc,
        ref_build=ref_build,
        ref_fasta=ref_fasta,
        bin_size=bin_size,
        athena_docker=athena_docker,
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_chrom_shard=runtime_attr_chrom_shard,
        runtime_attr_make_pairs=runtime_attr_make_pairs,
        runtime_attr_make_pairs=runtime_attr_annotate_pairs,
        runtime_attr_merge_annotated_pairs=runtime_attr_merge_annotated_pairs
    }

    # Step 4. Uniformly sample pairs per chromosome

    # (Step 5 called outside of scatter; see below)

    # Step 6. Apply PCA transformation to 2D pairs


    # Step 7. Merge & sort all 2D bins per chromosome

  }

  # [Optional] Step 5.1. Visualize distributions of all features
  # if visualize_features {
  #   call VisualizeFeatures {
  #     input:
  #   }
  # }

  # Step 5.2. Learn PCA transformation from sampled bins

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
      --exclusion-list-all ~{bin_exclusion_mask} \
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
    File bedtools_genome_file = "contigs.genome"
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

