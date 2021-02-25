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
    Int 
  }

}