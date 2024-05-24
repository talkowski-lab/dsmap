#######################
#    DSMap Project    #
#######################
#
# Utils.wdl
#
# Small utility tasks reused across workflows
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Structs.wdl"


# Subset a tabix-indexed input file to a single chromosome
# and divide that subset into smaller shard of prespecified size
task SingleChromShard {
  input {
    File infile
    File infile_idx
    String contig
    Int shard_size
    String prefix
    String file_format

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size(infile, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    # Extract header
    tabix -H ~{infile} > header

    # Extract chromosome of interest and split into smaller shards
    tabix ~{infile} "~{contig}" \
    | split -d -a 7 -l ~{shard_size} - shardfile_
    find ./ -name "shardfile_*"

    # Clean up shards
    n_shards=$( find ./ -name "shardfile_*" | wc -l )
    for i in $( seq -w 0000000 $(( $n_shards - 1 )) ); do
      cat header shardfile_$i \
      | bgzip -c \
      > "~{prefix}.~{contig}.shard_$i.~{file_format}.gz"
      # Tabix to ensure the shards are properly formatted
      tabix -p "~{file_format}" -f \
        "~{prefix}.~{contig}.shard_$i.~{file_format}.gz"
    done
  }

  output {
    Array[File] shards = glob("~{prefix}.~{contig}.shard_*.~{file_format}.gz")
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


# Merge, sort, and tabix an array of BED files
task MergeBEDs {
  input {
    Array[File] beds
    String prefix
    Boolean beds_are_bgzipped = true

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 10 + ceil(2 * size(beds, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    if [ ~{beds_are_bgzipped} == "true" ]; then
      zcat ~{beds[0]} | grep -e '^#' > header.tsv
      zcat ~{sep=" " beds} \
      | grep -ve '^#' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.tsv - \
      | bgzip -c \
      > "~{prefix}.bed.gz"
    else
      grep -e '^#' ~{beds[0]} > header.tsv
      cat ~{sep=" " beds} \
      | grep -ve '^#' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.tsv - \
      | bgzip -c \
      > "~{prefix}.bed.gz"
    fi

    tabix -p bed -f "~{prefix}.bed.gz"
  }

  output {
    File merged_bed = "~{prefix}.bed.gz"
    File merged_bed_idx = "~{prefix}.bed.gz.tbi"
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


# Visualize feature distributions from a BED3+
task VisualizeFeatures {
  input {
    File bed
    File? transformations_tsv
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 10 + ceil(2 * size(bed, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    mkdir ~{prefix}_feature_hists

    # Create option for feature transformation prior to model learning
    if [ "~{defined(transformations_tsv)}" == "true" ]; then
      athena_options="--transformations-tsv ~{transformations_tsv}"
    else
      athena_options=""
    fi

    # Visualize features with athena
    athena_cmd="athena feature-hists $athena_options --ignore-columns 3"
    athena_cmd="$athena_cmd ~{bed} ~{prefix}_feature_hists/~{prefix}"
    echo -e "Now visualizing using command:\n$athena_cmd"
    eval $athena_cmd
    tar -czvf "~{prefix}_feature_hists.tar.gz" "~{prefix}_feature_hists"
  }

  output {
    File feature_hists = "~{prefix}_feature_hists.tar.gz"
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


# Filter an input BED versus an exclusion bed
task ApplyExclusionBED {
  input {
    File inbed
    File exbed
    String prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10 + ceil(2 * size(inbed, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    bedtools intersect -v -header -wa \
      -a ~{inbed} \
      -b ~{exbed} \
    | bgzip -c \
    > ~{prefix}.bed.gz
    tabix -f ~{prefix}.bed.gz
  }

  output {
    File filtered_bed = "~{prefix}.bed.gz"
    File filtered_bed_idx = "~{prefix}.bed.gz.tbi"
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


# Compress an aribtrary set of files into a single tarball for easy I/O
task MakeTarball {
  input {
    Array[File] files_to_tar
    String tarball_prefix

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: 10 + ceil(2 * size(files_to_tar, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    mkdir ~{tarball_prefix}

    mv ~{sep=" " files_to_tar} ~{tarball_prefix}/

    tar -czvf ~{tarball_prefix}.tar.gz ~{tarball_prefix}
  }

  output {
    File tarball = "~{tarball_prefix}.tar.gz"
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


# Plots a distribution of CNV mutation rates over all targets/query entities
task PlotMuDist {
  input {
    File mu_tsv
    String cnv
    Boolean distance = false
    String? x_title
    String? y_title
    String out_prefix

    String dsmap_r_docker

    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10 + ceil(2 * size(mu_tsv, "GB")),
    boot_disk_gb: 15,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Build & execute command
    plot_cmd="/opt/dsmap/scripts/mu/plot_mu_distrib.R --cnv ~{cnv}"
    if [ ~{defined(x_title)} == "true" ]; then
      plot_cmd="$plot_cmd --x-title \"~{x_title}\""
    fi
    if [ ~{defined(y_title)} == "true" ]; then
      plot_cmd="$plot_cmd --y-title \"~{y_title}\""
    fi
    plot_cmd="$plot_cmd ~{true='--distance' false='' distance} ~{mu_tsv} ~{out_prefix}"
    echo -e "Now plotting mutation rate histogram using command:\n$plot_cmd"
    eval $plot_cmd
  >>>

  output {
    File mu_dist = "~{out_prefix}.mutation_rate_~{true='vs_distance' false='hist' distance}.pdf"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: dsmap_r_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
