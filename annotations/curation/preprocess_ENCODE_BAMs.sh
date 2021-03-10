#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to preprocess & curate all ENCODE BAMs


# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-db052a-f527e6

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Prep Cromwell directory structure
mkdir cromwell
for subdir in logs outputs inputs; do
  mkdir cromwell/$subdir
done

# Make inputs .json
cat <<EOF > cromwell/inputs/ENCODE.PreprocessExternalBAMs.input.json
{
  "PreprocessExternalBAMs.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:latest",
  "PreprocessExternalBAMs.bam_urls_tsv": "gs://dsmap/data/annotations/ENCODE/bams/ENCODE.bams_to_curate.tsv"
}
EOF

# Make dependencies .zip
if [ -e dsmap.dependencies.zip ]; then
  rm dsmap.dependencies.zip
fi
cd /opt/dsmap/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip / && \
cd -

# Submit to cromwell
cd /opt/dsmap && \
cromshell -t 20 submit \
  wdls/PreprocessExternalBAMs.wdl \
  /cromwell/inputs/ENCODE.PreprocessExternalBAMs.input.json \
  cromwell/dsmap.cromwell_options.json \
  /dsmap.dependencies.zip \
> /cromwell/logs/ENCODE.PreprocessExternalBAMs.log && \
cd -

# Check status
cromshell -t 20 metadata \
  $( cat /cromwell/logs/ENCODE.PreprocessExternalBAMs.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status
