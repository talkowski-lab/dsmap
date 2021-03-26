#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Development code for prototype WGS mutation rate model
# SV intersection and mutation rate model training & evaluation


###############
#    Setup    #
###############
# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-bb809b-33aeba

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap
gcloud auth application-default login

# Set global parameters
export cohort="HGSV"
export tech="WGS"
export cnv="DEL"
export vcf="gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.wCI.vcf.gz"
export vcf_idx="gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.wCI.vcf.gz.tbi"
export main_prefix="${cohort}.${tech}.dev"


############################################################
#    Build mutation rate model - step by step, locally     #
############################################################
# Note: in practice, this will be conducted in parallel for all contigs in master
# workflow, but for local dev purposes we are restricting to a single contig
export contig="chr22"

# Download necessary inputs to TrainMuModel.wdl
gsutil -m cp \
  gs://dsmap/data/dev/$main_prefix.BinAndAnnotateGenome/$main_prefix.pairs.eigen.$contig.bed.gz* \
  $vcf_idx \
  ./


################################################
#    Step 1. Intersect CNVs with bin pairs     #
################################################
# Note: in practice, this will be conducted in parallel for all contigs in master
# workflow, but for local dev purposes we are restricting to a single contig
export pairs_bed="$main_prefix.pairs.eigen.$contig.bed.gz"
export pairs_bed_idx="$main_prefix.pairs.eigen.$contig.bed.gz.tbi"
export prefix="$main_prefix.$cnv"

# Localize variants from contig
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
tabix -h $vcf $contig | bgzip -c > $prefix.$contig.svs.vcf.gz

# Intersect variants with pairs
athena count-sv \
  --bin-format 2D \
  --comparison breakpoint \
  --probabilities \
  --bgzip \
  $pairs_bed \
  $prefix.$contig.svs.vcf.gz \
  $prefix.pairs_wCounts.bed.gz
tabix -f $prefix.pairs_wCounts.bed.gz





