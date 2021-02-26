#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Development code for prototype WGS mutation rate model


###############
#    Setup    #
###############
# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:latest

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Set global parameters
export cohort="HGSV"
export tech="WGS"
export prefix="${cohort}.${tech}.dev"
export athena_docker="us.gcr.io/broad-dsmap/athena:latest"

# Prep Cromwell directory structure
mkdir cromwell
for subdir in logs outputs inputs; do
  mkdir cromwell/$subdir
done


#############################
#    Prep training data     #
#############################
# Download filtered HGSV VCFs
gsutil -m cp gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.*.sites.vcf.gz* ./

# Get CNV distributions
for CNV in DEL DUP; do
  echo -e "\n\n$CNV:"
  athena vcf-stats $cohort.$tech.filtered.$CNV.sites.vcf.gz
done


######################################################
#    Bin and annotate genome - one-shot workflow     #
######################################################
# Make inputs .json
cat <<EOF > cromwell/inputs/$prefix.BinAndAnnotateGenome.input.json
{
  "BinAndAnnotateGenome.bin_exclusion_mask": "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.athena_docker": "us.gcr.io/broad-dsmap/athena:wgs-mu-dev",
  "BinAndAnnotateGenome.bins_per_shard": 1000,
  "BinAndAnnotateGenome.prefix": "$prefix",
  "BinAndAnnotateGenome.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
  "BinAndAnnotateGenome.bin_size": 25000
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
  wdls/BinAndAnnotateGenome.wdl \
  /cromwell/inputs/$prefix.BinAndAnnotateGenome.input.json \
  cromwell/dsmap.cromwell_options.json \
  /dsmap.dependencies.zip \
> /cromwell/logs/$prefix.BinAndAnnotateGenome.log && \
cd -

# Check status
cromshell -t 20 metadata \
  $( cat /cromwell/logs/$prefix.BinAndAnnotateGenome.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status


##########################################################
#    Bin and annotate genome - step by step, locally     #
##########################################################
# Download necessary inputs to BinAndAnnotateGenome.wdl
gsutil -m cp \
  gs://dsmap/data/dev/hg38.contigs.dev.fai \
  gs://dsmap/data/references/hg38.nmask.bed.gz \
  ./

# Set parameters as required in BinAndAnnotateGenome.wdl
export contigs_fai=hg38.contigs.dev.fai
export bin_exclusion_mask=hg38.nmask.bed.gz
export bin_size=25000
export bins_per_shard=1000

# Step 1. Create all 1D bins
cut -f1-2 ${contigs_fai} > contigs.genome
athena make-bins \
  --exclusion-list-all ${bin_exclusion_mask} \
  --buffer ${bin_size} \
  --include-chroms $( cut -f1 contigs.genome | paste -s -d, ) \
  --bgzip \
  contigs.genome \
  ${bin_size} \
  ${prefix}.bins.bed.gz
tabix -f ${prefix}.bins.bed.gz

# Step 2: annotate 1D bins per contig
# (Note: for dev. purposes, subsetting this to 1,000 bins on chr22 as 
#  it would be handled in BinAndAnnotateGenome.wdl)
export contig=chr22
cat <( tabix -H ${prefix}.bins.bed.gz ) \
    <( tabix ${prefix}.bins.bed.gz chr22 | head -n${bins_per_shard} ) \
| bgzip -c \
> ${prefix}.${contig}.shard_000000.bed.gz









