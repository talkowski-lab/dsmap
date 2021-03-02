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

import "Utils.wdl"
import "Structs.wdl"

workflow BinAndAnnotateGenome {
	input {
    File contigs_fai
    # File pair_annotations_list
    File bin_exclusion_mask
    String prefix
    Int bin_size
    # Int max_pair_distance
    Int bins_per_shard
    # Int pairs_for_pca

    File? bin_annotations_list
    File? bin_annotations_list_remote_tabix
    File? bin_annotations_list_ucsc
    String? ref_build = "hg38"
    File? ref_fasta
    File? ref_fasta_idx

    String athena_docker
    String athena_cloud_docker

    RuntimeAttr? runtime_attr_make_bins
    RuntimeAttr? runtime_attr_chrom_shard_1d
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

    # Shard bins from each chromosome
    call Utils.SingleChromShard as ChromShard1D {
      input:
        infile=MakeBins.bins_bed,
        infile_idx=MakeBins.bins_bed_idx,
        contig=contig[0],
        shard_size=bins_per_shard,
        prefix=prefix,
        file_format="bed",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_chrom_shard_1d
    }

    # Scatter over shards per chromosome
    scatter( shard in ChromShard1D.shards ) {
    
      # Step 2. Annotate 1D bins per contig
      call AnnotateBins1D {
        input:
          bed=shard,
          bin_annotations_list=bin_annotations_list,
          bin_annotations_list_remote_tabix=bin_annotations_list_remote_tabix,
          bin_annotations_list_ucsc=bin_annotations_list_ucsc,
          ref_build=ref_build,
          ref_fasta=ref_fasta,
          ref_fasta_idx=ref_fasta_idx,
          athena_cloud_docker=athena_cloud_docker
      }

      # Step 3. Pair 2D bins and add 2D annotations

      # Step 4. Uniformly sample pairs per chromosome

      # (Step 5 called outside of scatter; see below)

      # Step 6. Apply PCA transformation to 2D pairs

    }

    # Step 7. Merge & sort all 2D bins per chromosome

  }

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


# Annotate a BED file of 1D bins
task AnnotateBins1D {
  input {
    File bed
    String ref_build
    
    File? bin_annotations_list
    File? bin_annotations_list_remote_tabix
    File? bin_annotations_list_ucsc
    File? ref_fasta
    File? ref_fasta_idx
    File? snv_mutrates_tsv

    String athena_cloud_docker
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

    # Prepare inputs
    out_prefix=$( echo "~{bed}" | sed 's/\.bed\.gz//g' )
    zcat ~{bed} | fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > regions.bed
    GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    # Slice large input files hosted remotely with athena slice-remote
    if [ ! -z ~{bin_annotations_list} ]; then
      cat ~{bin_annotations_list} > local_tracks.tsv
    fi
    if [ ! -z ~{bin_annotations_list_remote_tabix} ]; then
      if [ ! -z ${ref_fasta} ]; then
        remote_options="--ref-fasta ~{ref_fasta}"
      else
        remote_options=""
      fi
      athena slice-remote $remote_options \
        ~{bin_annotations_list_remote_tabix} \
        regions.bed \
      >> local_tracks.tsv
    fi

    # Build options for athena annotate-bins
    athena_options=""
    if [ ! -z local_tracks.tsv ]; then
      athena_options="$athena_options --track-list local_tracks.tsv"
    fi
    if [ ! -z ~{bin_annotations_list_ucsc} ]; then
      athena_options="$athena_options --ucsc-list ~{bin_annotations_list_ucsc}"
    fi
    if [ ! -z ~{ref_fasta} ]; then
      athena_options="$athena_options --fasta ~{ref_fasta}"
    fi
    if [ ! -z ~{snv_mutrates_tsv} ]; then
      athena_options="$athena_options --snv-mutrate ~{snv_mutrates_tsv}"
    fi

    # Annotate bins with athena
    athena_cmd="athena annotate-bins --ucsc-ref ~{ref_build} $athena_options"
    athena_cmd="$athena_cmd --no-ucsc-chromsplit --bgzip"
    athena_cmd="$athena_cmd ~{bed} ${out_prefix}.annotated.bed.gz"
    echo -e "Now annotating using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ${out_prefix}.annotated.bed.gz
  }

  output {
    File annotated_bed = "*.annotated.bed.gz"
    File annotated_bed_idx = "*.annotated.bed.gz.tbi"
  }
}
