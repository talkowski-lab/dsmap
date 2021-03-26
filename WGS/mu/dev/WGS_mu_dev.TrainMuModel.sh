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
export main_prefix="${cohort}.${tech}.dev"

# Prep Cromwell directory structure
mkdir cromwell
for subdir in logs outputs inputs; do
  mkdir cromwell/$subdir
done

########################################################
#    Train mutation rate model - one-shot workflow     #
########################################################
# Make dependencies .zip
if [ -e dsmap.dependencies.zip ]; then
  rm dsmap.dependencies.zip
fi
cd /opt/dsmap/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip / && \
cd -

# Must run once for each CNV type
for cnv in DEL DUP; do
  # Make inputs .json
  cat <<EOF > cromwell/inputs/$main_prefix.TrainMuModel.$cnv.input.json
{
  "TrainMuModel.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:latest",
  "TrainMuModel.cnv": "$cnv",
  "TrainMuModel.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
  "TrainMuModel.pairs_bed_prefix": "$main_prefix.pairs.eigen",
  "TrainMuModel.pairs_bucket": "gs://dsmap/data/dev/$main_prefix.BinAndAnnotateGenome",
  "TrainMuModel.prefix": "$main_prefix",
  "TrainMuModel.vcf": "gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.vcf.gz",
  "TrainMuModel.vcf_idx": "gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.vcf.gz"
}
EOF

  # Submit to cromwell
  cd /opt/dsmap && \
  cromshell -t 20 submit \
    wdls/TrainMuModel.wdl \
    /cromwell/inputs/$main_prefix.TrainMuModel.$cnv.input.json \
    cromwell/dsmap.cromwell_options.json \
    /dsmap.dependencies.zip \
  > /cromwell/logs/$main_prefix.TrainMuModel.$cnv.log && \
  cd -
done

# Check status
for cnv in DEL DUP; do
  cromshell -t 20 metadata \
    $( cat /cromwell/logs/$main_prefix.TrainMuModel.$cnv.log | tail -n1 | jq .id | tr -d '"' ) \
  | jq .status
done

# Sanity check: count the number of bins with SVs per chromosome
gsutil -m cp gs://dsmap/data/dev/hg38.contigs.dev.fai ./
for cnv in DEL DUP; do
  bucket=$( cromshell -t 20 metadata \
              $( cat /cromwell/logs/$main_prefix.TrainMuModel.$cnv.log | tail -n1 | jq .id | tr -d '"' ) \
            | jq .workflowRoot | tr -d '"' )
  while read shard contig; do
    echo -e "\n\n__$cnv on ${contig}__"
    gsutil -m cat ${bucket}call-IntersectSVs/shard-$shard/HGSV.WGS.dev.$cnv.pairs_wCounts.bed.gz \
    | gunzip -c | sed '1d' | cut -f4 | sort -n | uniq -c
  done < <( awk -v OFS="\t" '{ print NR-1, $1 }' hg38.contigs.dev.fai )
done

# TODO: Move all final outputs to dsmap gs:// bucket for permanent storage



############################################################
#    Train mutation rate model - step by step, locally     #
############################################################
# Note: in practice, this will be conducted in parallel for all contigs in master
# workflow, and will be launched once each for deletions and duplications, but 
# for local dev purposes we are restricting to deletions from a single contig
export contig="chr22"
export cnv="DEL"
export vcf="gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.vcf.gz"
export vcf_idx="gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.vcf.gz.tbi"

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
  $prefix.$pairs.wCounts.contig.bed.gz
tabix -f $prefix.$pairs.wCounts.contig.bed.gz





