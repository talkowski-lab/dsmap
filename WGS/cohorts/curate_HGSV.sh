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
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-7411a1-722506

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Set global parameters
export sv_pipeline_docker="us.gcr.io/broad-dsmap/gatk-sv/sv-pipeline:rlc_posthoc_filtering_cnv_mcnv_compatability_b874f7"
export sv_base_mini_docker="us.gcr.io/broad-dsmap/gatk-sv/sv-base-mini:rlc_posthoc_filtering_cnv_mcnv_compatability_b874f7"
export cohort="HGSV"
export tech="WGS"
export prefix="${cohort}.${tech}"
export raw_vcf="gs://dsmap/data/WGS/HGSV/unfiltered_full_vcf/HGSV.WGS.vcf.gz"

# Prep Cromwell directory structure
mkdir cromwell
for subdir in logs outputs inputs; do
  mkdir cromwell/$subdir
done

# Clone GATK-SV
git clone https://github.com/broadinstitute/gatk-sv.git /opt/gatk-sv


######################################################
#    Remove existing allele frequency annotations    #
######################################################
# Make inputs .json
cat <<EOF > cromwell/inputs/$prefix.CleanAFInfo.input.json
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
cromshell -t 20 submit \
  wdls/CleanAFInfo.wdl \
  /cromwell/inputs/$prefix.CleanAFInfo.input.json \
  cromwell/dsmap.cromwell_options.json \
  /dsmap.dependencies.zip \
> /cromwell/logs/$prefix.CleanAFInfo.log && \
cd -

# Check status
cromshell -t 20 metadata \
  $( cat /cromwell/logs/$prefix.CleanAFInfo.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status

# Write output paths to file
cromshell -t 20 metadata \
  $( cat /cromwell/logs/$prefix.CleanAFInfo.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .outputs | awk '{ print $2 }' | tr -d '"' | sed 's/,$//g' | sed '/^$/d' \
> /cromwell/outputs/$prefix.CleanAFInfo.outputs.txt


##########################################################################################
# Remove children from trios and re-annotate allele frequencies for all superpopulations #
##########################################################################################
# Collect necessary input files
vcf_noFreqs=$( grep -e '\.vcf.gz$' /cromwell/outputs/$prefix.CleanAFInfo.outputs.txt )
gsutil -m cp gs://dsmap/data/WGS/HGSV/$prefix.pop_labels.tsv ./
# Note: NA24385 is an isolated individual of Ashkenazi ancestry and also needs to be pruned from the VCF
gsutil -m cat gs://dsmap/data/WGS/HGSV/$prefix.ped \
| awk -v FS="\t" '{ if ($3!="0" || $4!="0") print $2 }' \
| cat - <( echo "NA24385" ) \
| sort -V | uniq > $prefix.samples_to_prune.txt
gsutil -m cp $prefix.samples_to_prune.txt gs://dsmap/data/WGS/HGSV/

# Make inputs .json
cat <<EOF > /cromwell/inputs/$prefix.PruneAndAddVafs.input.json
{
  "PruneAndAddVafs.vcf": "$vcf_noFreqs",
  "PruneAndAddVafs.sv_base_mini_docker": "$sv_base_mini_docker",
  "PruneAndAddVafs.ped_file": "gs://dsmap/data/WGS/HGSV/$prefix.ped",
  "PruneAndAddVafs.prune_list": "gs://dsmap/data/WGS/HGSV/$prefix.samples_to_prune.txt",
  "PruneAndAddVafs.sample_pop_assignments": "gs://dsmap/data/WGS/HGSV/$prefix.pop_labels.tsv",
  "PruneAndAddVafs.sv_pipeline_docker": "$sv_pipeline_docker",
  "PruneAndAddVafs.sv_per_shard": 5000,
  "PruneAndAddVafs.contig_list": "gs://dsmap/data/references/hg38.contigs.fai",
  "PruneAndAddVafs.prefix": "$prefix",
  "PruneAndAddVafs.vcf_idx": "${vcf_noFreqs}.tbi"
}
EOF

# Make dependencies .zip
if [ -e gatksv.dependencies.zip ]; then
  rm gatksv.dependencies.zip
fi
cd /opt/gatk-sv/wdl/ && \
zip gatksv.dependencies.zip *wdl && \
mv gatksv.dependencies.zip / && \
cd -

# Submit to cromwell
cd /opt/gatk-sv && \
cromshell -t 20 submit \
  wdl/PruneAndAddVafs.wdl \
  /cromwell/inputs/$prefix.PruneAndAddVafs.input.json \
  /opt/dsmap/cromwell/dsmap.cromwell_options.no_caching.json \
  /gatksv.dependencies.zip \
> /cromwell/logs/$prefix.PruneAndAddVafs.log && \
cd -

# Check status
cromshell -t 20 metadata \
  $( cat /cromwell/logs/$prefix.PruneAndAddVafs.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status

# Write output paths to file
cromshell -t 20 metadata \
  $( cat /cromwell/logs/$prefix.PruneAndAddVafs.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .outputs | awk '{ print $2 }' | tr -d '"' | sed 's/,$//g' | sed '/^$/d' \
> /cromwell/outputs/$prefix.PruneAndAddVafs.outputs.txt


###########################################################
# Filter to high-quality, rare, biallelic, autosomal CNVs #
###########################################################
# Download the VCF with new AF annotations and cut to sites.vcf only
vcf_wFreqs=$( grep -e '\.vcf.gz$' /cromwell/outputs/$prefix.PruneAndAddVafs.outputs.txt )
gsutil -m cat $vcf_wFreqs \
| gunzip -c | cut -f1-8 | bgzip -c \
> $prefix.wAFs.sites.vcf.gz
tabix -f $prefix.wAFs.sites.vcf.gz

# Save a copy of the annotated VCF for archival purposes in main DSMap storage bucket
gsutil -m cp \
  $prefix.wAFs.sites.vcf.gz \
  $prefix.wAFs.sites.vcf.gz.tbi \
  gs://dsmap/data/WGS/HGSV/unfiltered_sites_vcf/

# Filter VCF
for CNV in DEL DUP; do
  athena vcf-filter \
    --include-chroms "$( seq 1 22 | awk '{ print "chr"$1 }' | paste -s -d, )" \
    --svtypes $CNV \
    --maxAF 0.01 \
    --af-field POPMAX_AF \
    --minAC 1 \
    --minAN 2331 \
    --minQUAL 2 \
    --pHWE 0.01 \
    --bgzip \
    $prefix.wAFs.sites.vcf.gz \
    $prefix.filtered.$CNV.sites.vcf.gz
  tabix -f $prefix.filtered.$CNV.sites.vcf.gz
done

# Copy filtered data to DSMap bucket
gsutil -m cp \
  $prefix.filtered.*.sites.vcf.gz* \
  gs://dsmap/data/WGS/HGSV/


##############################
# Add breakpoint uncertainty #
##############################
# Add breakpoint uncertainty to VCFs
# Note: for development purposes, just applying a universal breakpoint uncertainty of ±50bp
# This will be revisited after developing GATK-SV uncertainty model (TODO)
for CNV in DEL DUP; do
  athena breakpoint-confidence \
    --min-ci 100 \
    --bgzip \
  $prefix.filtered.$CNV.sites.vcf.gz \
  $prefix.filtered.$CNV.sites.wCI.vcf.gz
  tabix -f $prefix.filtered.$CNV.sites.wCI.vcf.gz
done

# Copy data with breakpoint confidence to DSMap bucket
gsutil -m cp \
  $prefix.filtered.*.sites.wCI.vcf.gz* \
  gs://dsmap/data/WGS/HGSV/

