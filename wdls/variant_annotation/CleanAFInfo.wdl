version 1.0

import "../Structs.wdl"

workflow CleanAFInfo {
  input {
    File vcf
    File vcf_idx
    String prefix
    String dsmap_docker
  }

  call cleanVCF {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      prefix=prefix,
      dsmap_docker=dsmap_docker
  }

  output {
    File cleaned_vcf = cleanVCF.cleaned_vcf
    File cleaned_vcf_idx = cleanVCF.cleaned_vcf_idx
  }
}


task cleanVCF {
  input {
    File vcf
    File vcf_idx
    String prefix
    String dsmap_docker
    RuntimeAttr? runtime_attr_override
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
    /opt/dsmap/scripts/variant_annotation/clean_af_info.py \
      ~{vcf} \
      stdout \
    | bgzip -c \
    > ~{prefix}.noAFs.vcf.gz
    tabix -f ${prefix}.noAFs.vcf.gz
  }

  output {
    File cleaned_vcf = "~{prefix}.noAFs.vcf.gz"
    File cleaned_vcf_idx = "~{prefix}.noAFs.vcf.gz.tbi"
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
