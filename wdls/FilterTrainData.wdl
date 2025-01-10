#######################
#    DSMap Project    #
#######################
#
# FilterTrainData.wdl
#
# Filter bin-pair annotations and counts to data to be used for model training
#
# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>


version 1.0

import "Utils.wdl"
import "Structs.wdl"


workflow FilterTrainData {
  input {
    # General inputs
    String annotated_pairs_w_counts_bucket
    String annotated_pairs_w_counts_bed_prefix
    Boolean counts_are_probs
    File contigs_fai
    File training_mask
    String cnv
    String prefix

    # Diagnostics options
    Boolean run_diagnostics = true

    # Dockers
    String athena_docker
    String dsmap_r_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_apply_training_mask
    RuntimeAttr? runtime_attr_diagnostics
  }

  Array[String] contigs = transpose(read_tsv(contigs_fai))[0]

  # Parallelize per chromosome
  scatter ( contig in contigs ) {
    # Infer pairs BED file path
    File annotated_pairs_w_counts_bed = (
        annotated_pairs_w_counts_bucket + "/" + annotated_pairs_w_counts_bed_prefix +
        "." + contig + ".bed.gz"
    )
    File annotated_pairs_w_counts_bed_idx = annotated_pairs_w_counts_bed + ".tbi"

    # Step 1. Apply training mask
    call Utils.ApplyExclusionBED as ApplyTrainingMask {
      input:
        inbed=annotated_pairs_w_counts_bed,
        exbed=training_mask,
        prefix=basename(annotated_pairs_w_counts_bed, ".bed.gz") + ".training",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_apply_training_mask
    }
  }

  # [Optional] Step 2. Run diagnostics
  if ( run_diagnostics ) {
    call Utils.GetPairDiagnostics {
        input:
        pair_counts=ApplyTrainingMask.filtered_bed,
        counts_are_probs=counts_are_probs,
        cnv=cnv,
        prefix=prefix,
        dsmap_r_docker=dsmap_r_docker,
        runtime_attr_override=runtime_attr_diagnostics
    }

    # Tar diagnostics for convenience
    call Utils.MakeTarball as MergeDiagnostics {
        input:
        files_to_tar=GetPairDiagnostics.outputs,
        tarball_prefix="~{prefix}.~{cnv}.FilterTrainData.diagnostics",
        athena_docker=athena_docker,
        runtime_attr_override=runtime_attr_diagnostics
    }
  }

  output {
    Array[File] training_beds = ApplyTrainingMask.filtered_bed
    Array[File] training_bed_idxs = ApplyTrainingMask.filtered_bed_idx
    File? diagnostics = MergeDiagnostics.tarball
  }
}
