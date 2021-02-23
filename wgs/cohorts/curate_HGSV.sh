#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to curate HGSV WGS callset for DSMap modeling


###############
#    Setup    #
###############
# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap:hgsv-wgs-curation

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Set global parameters
export sv_pipeline_docker="us.gcr.io/talkowski-sv-gnomad/rlc/sv-pipeline:rlc_posthoc_filtering_cnv_mcnv_compatability_fb8bc8"
export cohort="HGSV"
export tech="WGS"
export prefix="${cohort}.${tech}"
export raw_vcf="gs://talkowski-sv-gnomad-output/1KGP/final_vcf/1KGP_2504_and_698_with_GIAB.boost.PASS_gt_revised.vcf.gz"

# Clone GATK-SV
git clone git@github.com:broadinstitute/gatk-sv.git /opt/gatk-sv


######################################################
#    Remove existing allele frequency annotations    #
######################################################
gsutil -m cp ${raw_vcf} ${raw_vcf}.tbi ./
/opt/dsmap/scripts/variant_annotation/clean_af_info.py \
  $( basename $raw_vcf ) \
  stdout \
| bgzip -c \
> $cohort.$tech.noAFs.vcf.gz
tabix -f $cohort.$tech.noAFs.vcf.gz
gsutil -m cp $cohort.$tech.noAFs.vcf.gz gs://dsmap/scratch/


##########################################################################################
# Remove children from trios and re-annotate allele frequencies for all superpopulations #
##########################################################################################


##################################################################
# Filter to high-quality, rare, biallelic, autosomal CNVs in HWE #
##################################################################
for CNV in DEL DUP; do
  athena vcf-filter -z \
    --exclude-chroms X,Y,M \
    --svtypes DEL \
    --maxAF 0.01 \
    --af-field POPMAX_AF \
    --minAC 1 \
    # --minAN 20402 \
    --minQUAL 100 \
    --pHWE 0.01 \
    gnomad_v2.1_sv.sites.vcf.gz \
    athena_training_deletions.vcf.gz

