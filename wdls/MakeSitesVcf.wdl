#######################
#    DSMap Project    #
#######################
#
# MakeSitesVcf.wdl
#
# Helper workflow to make a sites VCF from one or more VCFs
#
# Copyright (c) 2023-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rcollins@broadinstitute.org>


version 1.0

import "Structs.wdl"


workflow MakeSitesVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    String athena_docker

    RuntimeAttr? runtime_attr_override
  }

  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {

    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right
    String out_filename = basename(vcf, ".vcf.gz") + ".sites.vcf.gz"

    call MakeSites {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        out_filename=out_filename,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_override
    }
  }

  output {
    Array[File] sites_vcfs = MakeSites.vcf_out
    Array[File] sites_vcf_idxs = MakeSites.vcf_idx_out
  }
}


task MakeSites {
  input {
    File vcf
    File vcf_idx
    String out_filename
    String athena_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: ceil(3 * size([vcf], "GB")) + 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {

    set -euo pipefail

    zcat ~{vcf} | cut -f1-8 | bgzip -c > ~{out_filename}
    tabix -p vcf -f ~{out_filename}

  }

  output {
    File vcf_out = "~{out_filename}"
    File vcf_idx_out = "~{out_filename}.tbi"
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
