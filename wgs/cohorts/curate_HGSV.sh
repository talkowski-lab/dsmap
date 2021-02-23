#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to curate HGSV WGS callset for DSMap modeling


# Run in Athena docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap


# Set global parameters
export sv_pipeline_docker="us.gcr.io/talkowski-sv-gnomad/rlc/sv-pipeline:rlc_posthoc_filtering_cnv_mcnv_compatability_fb8bc8"
export cohort="HGSV"
export tech="WGS"
export prefix="${cohort}.${tech}"
export raw_vcf="gs://talkowski-sv-gnomad-output/1KGP/final_vcf/1KGP_2504_and_698_with_GIAB.boost.PASS_gt_revised.vcf.gz"


# Remove existing allele frequency annotations
gsutil -m cp $raw_vcf ./
zcat $( basename $raw_vcf ) \
| fgrep -v "" \
| fgrep -v "" \
| 


# Remove children from trios and re-annotate allele frequencies for all superpopulations


# Filter to high-quality, rare, biallelic, autosomal CNVs in HWE

