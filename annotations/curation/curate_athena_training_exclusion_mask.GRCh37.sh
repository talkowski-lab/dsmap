#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to curate athena neutral sites training filter for GRCh37


# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-4c74d9-097434

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Download gnomAD transcript constraint data
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz

# Download Gencode v19 GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# Download GRCh37 GERP elements from Ensembl v75 (last GRCh37 release)
wget -O - http://ftp.ensembl.org/pub/release-75/bed/37way_eutherian_mammals.homo_sapiens.bed.gz \
| gunzip -c | sed '1d' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | bgzip -c \
> gerp_constrained_elements.homo_sapiens.bed.gz
tabix -f gerp_constrained_elements.homo_sapiens.bed.gz

# Create athena neutral sites mask
/opt/dsmap/scripts/annotations/make_athena_training_mask.py \
  --gtf gencode.v19.annotation.gtf.gz \
  --constraint gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
  --gerp gerp_constrained_elements.homo_sapiens.bed.gz \
  --outfile athena.GRCh37.training_exclusion.bed
sed 's/^chr//g' athena.GRCh37.training_exclusion.bed | bgzip -c \
> athena.GRCh37.training_exclusion.bed.gz

# Copy neutral sites mask to Google Bucket
gsutil -m cp \
  athena.GRCh37.training_exclusion.bed.gz \
  gs://dsmap/data/athena/
