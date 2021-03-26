#######################
#    DSMap Project    #
#######################
#
# TrainMuModel.wdl
#
# Intersect SVs with bin-pairs and train mutation rate model (with optional evaluation)
#
# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow TrainMuModel {
  input {
    # General inputs
    String vcf
    File vcf_idx
    String pairs_bucket
    String? pairs_bed_prefix
    File contigs_fai
    String prefix
    String cnv

    # Dockers
    String athena_cloud_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_intersect_svs
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
  scatter ( contig in contigs ) {

    # Infer pairs BED file path
    String pairs_filename = if defined(pairs_bed_prefix) then pairs_bed_prefix else prefix + ".pairs.eigen"
    File pairs_bed = pairs_bucket + "/" + pairs_filename + "." + contig + ".bed.gz" 
    File pairs_bed_idx = pairs_bed + ".bed.gz.tbi" 

    # Step 1. Intersect CNVs with bin-pairs
    call IntersectSVs {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        pairs_bed=pairs_bed,
        pairs_bed_idx=pairs_bed_idx,
        contig=contig,
        prefix="~{prefix}.~{cnv}",
        athena_cloud_docker=athena_cloud_docker,
        runtime_attr_override=runtime_attr_intersect_svs
    }
  }
}


task IntersectSVs {
  input {
    String vcf #VCF passed as string to allow for remote tabixing without localizing entire VCF
    File vcf_idx
    File pairs_bed
    File pairs_bed_idx
    String contig
    String prefix
    
    String athena_cloud_docker

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

    # Localize variants from contig
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    tabix -h $vcf $contig | bgzip -c > $prefix.$contig.vcf.gz

    

  }

  output {}
  
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
