#######################
#    DSMap Project    #
#######################
#
# FilterMu.wdl
#
# Filter mutation matrix
#
# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow FilterMu {
  input {
    # General inputs
    String mu_bucket
    String mu_bed_prefix
    String pairs_bucket
    String pairs_bed_prefix
    File contigs_fai
    String filter_prefix

    # Dockers
    String athena_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_filter_mu
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
  scatter ( contig in contigs ) {
    String del_prefix = ".DEL" + "." + contig
    String dup_prefix = ".DUP" + "." + contig

    File pairs_bed = pairs_bucket + "/" + pairs_bed_prefix + "." + contig + ".bed.gz"
    File del_mu_bed = mu_bucket + "/" + mu_bed_prefix + del_prefix + ".mu.bed.gz"
    File dup_mu_bed = mu_bucket + "/" + mu_bed_prefix + dup_prefix + ".mu.bed.gz"

    # Filter DEL and DUP mutation matrices to input pairs
    call Utils.ApplyMatchBED as FilterDelMuToPairs {
      input:
        inbed=del_mu_bed,
        matchbed=pairs_bed,
        prefix=mu_bed_prefix + "." + filter_prefix + del_prefix + ".mu.bed.gz",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_filter_mu
    }
    call Utils.ApplyMatchBED as FilterDupMuToPairs {
      input:
        inbed=dup_mu_bed,
        matchbed=pairs_bed,
        prefix=mu_bed_prefix + "." + filter_prefix + dup_prefix + ".mu.bed.gz",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_filter_mu
    }
  }

  output {
    Array[File] filtered_del_mu_beds = FilterDelMuToPairs.filtered_bed
    Array[File] filtered_del_mu_bed_idxs = FilterDelMuToPairs.filtered_bed_idx
    Array[File] filtered_dup_mu_beds = FilterDupMuToPairs.filtered_bed
    Array[File] filtered_dup_mu_bed_idxs = FilterDupMuToPairs.filtered_bed_idx
  }
}
