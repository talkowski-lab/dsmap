#######################
#    DSMap Project    #
#######################
#
# FilterVcf.wdl
#
# Filter VCF to high-quality rare CNVs
#
# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>


version 1.0

import "Structs.wdl"


workflow FilterVcf {
  input {
    # General inputs
    File vcf
    File vcf_idx
    Float? max_af
    Array[String] af_fields = []
    Int? min_ac
    Int? min_an
    Int? min_qual
    Float? min_p_hwe
    Array[File]? exclusion_beds
    Boolean split_cnv = true
    String? filter_prefix

    # Diagnostics options
    Boolean run_diagnostics = true

    # Dockers
    String athena_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf, ".vcf.gz") + if defined(filter_prefix) then "." + filter_prefix else ".filtered"

  if ( split_cnv ) {
    # Filter VCF to high-quality rare DELs and DUPs
    call FilterCnvs as FilterDels {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        cnv="DEL",
        max_af=max_af,
        af_fields=af_fields,
        min_ac=min_ac,
        min_an=min_an,
        min_qual=min_qual,
        p_hwe=min_p_hwe,
        exclusion_beds=exclusion_beds,
        prefix=prefix + ".DEL",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_override
    }
    call FilterCnvs as FilterDups {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        cnv="DUP",
        max_af=max_af,
        af_fields=af_fields,
        min_ac=min_ac,
        min_an=min_an,
        min_qual=min_qual,
        p_hwe=min_p_hwe,
        exclusion_beds=exclusion_beds,
        prefix=prefix + ".DUP",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_override
    }

    if ( run_diagnostics ) {
      # Compute stats on CNV size and spacing in DEL and DUP VCFs
      call GetVcfStats as GetDelStats {
        input:
          vcf=FilterDels.vcf_out,
          vcf_idx=FilterDels.vcf_idx_out,
          prefix=prefix + ".DEL",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_override
      }
      call GetVcfStats as GetDupStats {
        input:
          vcf=FilterDups.vcf_out,
          vcf_idx=FilterDups.vcf_idx_out,
          prefix=prefix + ".DUP",
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_override
      }
    }
  }
  
  if ( !split_cnv ) {
    # Filter VCF to high-quality rare CNVs
    call FilterCnvs {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        max_af=max_af,
        af_fields=af_fields,
        min_ac=min_ac,
        min_an=min_an,
        min_qual=min_qual,
        p_hwe=min_p_hwe,
        exclusion_beds=exclusion_beds,
        prefix=prefix,
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_override
    }

    if ( run_diagnostics ) {
      # Compute stats on CNV size and spacing in VCF
      call GetVcfStats {
        input:
          vcf=FilterCnvs.vcf_out,
          vcf_idx=FilterCnvs.vcf_idx_out,
          prefix=prefix,
          athena_docker=athena_docker,
          runtime_attr_override=runtime_attr_override
      }
    }
  }

  output {
    # Outputs if split_cnv is optioned
    File? del_vcf = FilterDels.vcf_out
    File? del_vcf_idx = FilterDels.vcf_idx_out
    File? dup_vcf = FilterDups.vcf_out
    File? dup_vcf_idx = FilterDups.vcf_idx_out
    File? del_vcf_stats = GetDelStats.stats_txt
    File? dup_vcf_stats = GetDupStats.stats_txt

    # Outputs if split_cnv is not optioned
    File? cnv_vcf = FilterCnvs.vcf_out
    File? cnv_vcf_idx = FilterCnvs.vcf_idx_out
    File? cnv_vcf_stats = GetVcfStats.stats_txt
  }
}


task FilterCnvs {
  input {
    File vcf
    File vcf_idx
    String? cnv
    String prefix

    Float? max_af
    Array[String] af_fields = []
    Int? min_ac
    Int? min_an
    Int? min_qual
    Float? p_hwe
    Array[File]? exclusion_beds

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

  command <<<

    set -euo pipefail

    # Create options for filtering
    athena_options=""
    if [ "~{defined(cnv)}" == "true" ]; then
      athena_options="$athena_options --svtypes ~{cnv}"
    fi
    if [ "~{defined(max_af)}" == "true" ]; then
      athena_options="$athena_options --maxAF ~{max_af}"
    fi
    if [ ~{length(af_fields)} -gt 0 ]; then
      athena_options="$athena_options --af-field ~{sep=' --af-field ' af_fields}"
    fi
    if [ "~{defined(min_ac)}" == "true" ]; then
      athena_options="$athena_options --minAC ~{min_ac}"
    fi
    if [ "~{defined(min_an)}" == "true" ]; then
      athena_options="$athena_options --minAN ~{min_an}"
    fi
    if [ "~{defined(min_qual)}" == "true" ]; then
      athena_options="$athena_options --minQUAL ~{min_qual}"
    fi
    if [ "~{defined(p_hwe)}" == "true" ]; then
      athena_options="$athena_options --pHWE ~{p_hwe}"
    fi
    if [ "~{defined(exclusion_beds)}" == "true" ]; then
      athena_options="$athena_options --exclusion-list ~{sep=' --exclusion-list ' exclusion_beds}"
    fi

    # Filter VCF with athena
    athena_cmd="athena vcf-filter $athena_options"
    athena_cmd="$athena_cmd --include-chroms '$( seq 1 22 | awk '{ print "chr"$0 }' | paste -s -d, )'"
    athena_cmd="$athena_cmd --bgzip ~{vcf} ~{prefix}.vcf.gz"
    echo -e "Now filtering VCF using command:\n$athena_cmd"
    eval $athena_cmd
    tabix -f ~{prefix}.vcf.gz

  >>>

  output {
    File vcf_out = "~{prefix}.vcf.gz"
    File vcf_idx_out = "~{prefix}.vcf.gz.tbi"
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


task GetVcfStats {
  input {
    File vcf
    File vcf_idx
    String prefix

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

    athena vcf-stats ~{vcf} -q '0.25,0.5,0.75,0.9,0.95,0.99,0.999,1.0' > ~{prefix}.stats.txt

  }

  output {
    File stats_txt = "~{prefix}.stats.txt"
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
