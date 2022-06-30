#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to curate athena neutral sites training filter for hg38


# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-7411a1-722506

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Download gnomAD transcript constraint data
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz

# Download Gencode v37 GTF
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz

# Download hg38 GERP elements from Ensembl v103
wget http://ftp.ensembl.org/pub/release-103/bed/ensembl-compara/90_mammals.gerp_constrained_element/gerp_constrained_elements.homo_sapiens.bb

# Convert GERP from BigBed to BED
cd opt/bin && \
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed && \
chmod a+x bigBedToBed && \
bigBedToBed \
  gerp_constrained_elements.homo_sapiens.bb \
  gerp_constrained_elements.homo_sapiens.bed
sort -Vk1,1 -k2,2n -k3,3n \
  gerp_constrained_elements.homo_sapiens.bed \
| awk '{ print "chr"$0 }' \
| bgzip -c \
> gerp_constrained_elements.homo_sapiens.bed.gz
tabix -f gerp_constrained_elements.homo_sapiens.bed.gz

# Create athena neutral sites mask
/opt/dsmap/scripts/annotations/make_athena_training_mask.py \
  --gtf gencode.v37.annotation.gtf.gz \
  --constraint gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
  --gerp gerp_constrained_elements.homo_sapiens.bed.gz \
  --outfile athena.hg38.training_exclusion.bed.gz \
  --bgzip

# Copy neutral sites mask to Google Bucket
gsutil -m cp \
  athena.hg38.training_exclusion.bed.gz \
  gs://dsmap/data/athena/
