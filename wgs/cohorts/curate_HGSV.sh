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
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:hgsv-wgs-curation

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Set global parameters
export sv_pipeline_docker="us.gcr.io/talkowski-sv-gnomad/rlc/sv-pipeline:rlc_posthoc_filtering_cnv_mcnv_compatability_fb8bc8"
export cohort="HGSV"
export tech="WGS"
export prefix="${cohort}.${tech}"
export raw_vcf="gs://talkowski-sv-gnomad-output/1KGP/final_vcf/1KGP_2504_and_698_with_GIAB.boost.PASS_gt_revised.vcf.gz"
export BASE_PATH="$PATH"
export CROM_CURL_PATH="/usr/bin:$PATH"

# Clone GATK-SV
git clone https://github.com/broadinstitute/gatk-sv.git /opt/gatk-sv


######################################################
#    Remove existing allele frequency annotations    #
######################################################
# Make inputs .json
cat <<EOF > $cohort.$tech.clean_af_info.input.json
{
  "CleanAFInfo.dsmap_docker": "us.gcr.io/broad-dsmap/dsmap:hgsv-wgs-curation",
  "CleanAFInfo.vcf_idx": "${raw_vcf}.tbi",
  "CleanAFInfo.vcf": "$raw_vcf",
  "CleanAFInfo.prefix": "$prefix"
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
export PATH=$CROM_CURL_PATH && \
cromshell -t 20 submit \
  wdls/CleanAFInfo.wdl \
  /$cohort.$tech.clean_af_info.input.json \
  cromwell/dsmap.cromwell_options.json \
  /dsmap.dependencies.zip \
> /CleanAFInfo.$cohort.$tech.cromshell_submission.log && \
export PATH=$BASE_PATH && \
cd -



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

