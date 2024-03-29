#######################
#    DSMap Project    #
#######################
#
# CleanAFInfo.wdl
#
# Pre-process remote hosted BAM files for athena annotation
#
# Input: three-column .tsv specifying BAM URL, desired filename prefix, and curated destination bucket
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Structs.wdl"

workflow PreprocessExternalBAMs {
  input {
    File bam_urls_tsv
    String athena_cloud_docker
    Boolean? check_for_prior_curation = true

    RuntimeAttr? runtime_attr
  }

  Array[Array[String]] bam_info = read_tsv(bam_urls_tsv)

  scatter ( bam_line in bam_info ) {
    call CurateBam {
      input:
        bam=bam_line[0],
        prefix=bam_line[1],
        dest_bucket=bam_line[2],
        check_for_prior_curation=check_for_prior_curation,
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_override=runtime_attr
    }
  }
}


# Download, sort, index, and upload BAM to specified destination bucket
task CurateBam {
  input {
    String bam
    String prefix
    String dest_bucket
    String athena_cloud_docker
    Boolean? check_for_prior_curation = true
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 4, 
    mem_gb: 16,
    disk_gb: 250,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    index_exists=$( gsutil ls ~{dest_bucket}/~{prefix}.bam.bai | wc -l || true )
    
    if [ "~{check_for_prior_curation}" == "true" ] && [ "$index_exists" == "1" ]; then
    
      echo "SKIPPING: found that BAM index already exists at ~{dest_bucket}/~{prefix}.bam.bai" > ~{prefix}.success.txt

    else

      wget -nv -O - ~{bam} \
      | samtools sort -@ 4 -m 2G -T ~{prefix}.tmpsort -O BAM -o ~{prefix}.bam -

      samtools index -b ~{prefix}.bam

      gsutil -m cp ~{prefix}.bam* ~{dest_bucket}/

      echo "SUCCESS: uploaded to ~{dest_bucket}/~{prefix}.bam" > ~{prefix}.success.txt

    fi
  }

  output {
    File success_token = "~{prefix}.success.txt"
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
