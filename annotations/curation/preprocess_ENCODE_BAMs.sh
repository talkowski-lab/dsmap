#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to preprocess & curate all ENCODE BAMs


# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:latest

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Download input file
gsutil cp gs://dsmap/data/annotations/ENCODE/bams/ENCODE.bams_to_curate.tsv ./

# # Make inputs .json
# cat <<EOF > cromwell/inputs/$prefix.BinAndAnnotateGenome.input.json
# {
#   "BinAndAnnotateGenome.bin_exclusion_mask": "gs://dsmap/data/references/hg38.nmask.bed.gz",
#   "BinAndAnnotateGenome.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
#   "BinAndAnnotateGenome.prefix": "$prefix",
#   "BinAndAnnotateGenome.bin_size": 25000,
#   "BinAndAnnotateGenome.bins_per_shard": 1000,
#   "BinAndAnnotateGenome.athena_docker": "us.gcr.io/broad-dsmap/athena:wgs-mu-dev",
#   "BinAndAnnotateGenome.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:wgs-mu-dev"
# }
# EOF

# # Make dependencies .zip
# if [ -e dsmap.dependencies.zip ]; then
#   rm dsmap.dependencies.zip
# fi
# cd /opt/dsmap/wdls/ && \
# zip dsmap.dependencies.zip *wdl && \
# mv dsmap.dependencies.zip / && \
# cd -

# # Submit to cromwell
# cd /opt/dsmap && \
# cromshell -t 20 submit \
#   wdls/BinAndAnnotateGenome.wdl \
#   /cromwell/inputs/$prefix.BinAndAnnotateGenome.input.json \
#   cromwell/dsmap.cromwell_options.json \
#   /dsmap.dependencies.zip \
# > /cromwell/logs/$prefix.BinAndAnnotateGenome.log && \
# cd -

# # Check status
# cromshell -t 20 metadata \
#   $( cat /cromwell/logs/$prefix.BinAndAnnotateGenome.log | tail -n1 | jq .id | tr -d '"' ) \
# | jq .status


