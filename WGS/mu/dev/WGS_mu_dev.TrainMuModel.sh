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
export docker_tag="wgs-mu-dev-58377c-ec4099"
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:$docker_tag

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap
gcloud auth application-default login

# Set global parameters
export docker_tag="wgs-mu-dev-58377c-ec4099"
export cohort="HGSV"
export tech="WGS"
export main_prefix="${cohort}.${tech}.dev"

# Prep Cromwell directory structure
mkdir cromwell
for subdir in logs outputs inputs; do
  mkdir cromwell/$subdir
done


#######################################################
#    Train mutation rate model - one-shot workflow    #
#######################################################
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
  cat <<EOF > /cromwell/inputs/$main_prefix.TrainMuModel.$cnv.input.json
{
  "TrainMuModel.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:$docker_tag",
  "TrainMuModel.athena_docker": "us.gcr.io/broad-dsmap/athena:$docker_tag",
  "TrainMuModel.cnv": "$cnv",
  "TrainMuModel.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
  "TrainMuModel.dsmap_r_docker" : "us.gcr.io/broad-dsmap/dsmap-r:$docker_tag",
  "TrainMuModel.model": "logit",
  "TrainMuModel.pairs_bed_prefix": "$main_prefix.pairs.eigen",
  "TrainMuModel.pairs_bucket": "gs://dsmap/data/dev/$main_prefix.BinAndAnnotateGenome",
  "TrainMuModel.prefix": "$main_prefix",
  "TrainMuModel.run_diagnostics" : true,
  "TrainMuModel.shard_size_apply_mu": 100000,
  "TrainMuModel.training_mask": "gs://dsmap/data/athena/athena.hg38.training_exclusion.bed.gz",
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

# Move all final outputs to dsmap gs:// bucket for permanent storage
for cnv in DEL DUP; do
  workflow_id=$( cat /cromwell/logs/$main_prefix.TrainMuModel.$cnv.log | tail -n1 | jq .id | tr -d '"' )
  cromshell -t 20 metadata $workflow_id \
  | jq .outputs | sed 's/\"/\n/g' | fgrep "gs://" \
  | gsutil -m cp -I gs://dsmap/data/dev/$main_prefix.TrainMuModel/
  # Also make note of most recent workflow ID (in order to reference intermediate 
  # outputs after Docker container is closed)
  echo "$workflow_id" > $main_prefix.$cnv.TrainMuModel.latest_workflow_id.txt
  gsutil -m cp \
    $main_prefix.$cnv.TrainMuModel.latest_workflow_id.txt \
    gs://dsmap/data/dev/$main_prefix.TrainMuModel/
done


###########################################################
#    Train mutation rate model - step by step, locally    #
###########################################################
# Note: in practice, this will be conducted in parallel for all contigs in master
# workflow, and will be launched once each for deletions and duplications, but 
# for local dev purposes we are restricting to deletions from a single contig
export contig="chr22"
export cnv="DEL"
export vcf="gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.vcf.gz"
export vcf_idx="gs://dsmap/data/$tech/$cohort/$cohort.$tech.filtered.$cnv.sites.vcf.gz.tbi"
export training_mask="athena.hg38.training_exclusion.bed.gz"

# Download necessary inputs to TrainMuModel.wdl
gsutil -m cp \
  gs://dsmap/data/dev/$main_prefix.BinAndAnnotateGenome/$main_prefix.pairs.eigen.$contig.bed.gz* \
  $vcf_idx \
  gs://dsmap/data/athena/athena.hg38.training_exclusion.bed.gz \
  gs://dsmap/data/dev/hg38.contigs.dev.fai \
  ./


###############################################
#    Step 1. Intersect CNVs with bin pairs    #
###############################################
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
  $prefix.pairs.wCounts.$contig.bed.gz
tabix -f $prefix.pairs.wCounts.$contig.bed.gz


#####################################
#    Step 2. Apply training mask    #
#####################################
# Set parameters as required for ApplyExclusionBED task in WDL
export inbed=$prefix.pairs.wCounts.$contig.bed.gz
export exbed=$training_mask
export prefix="$prefix.pairs.wCounts.$contig.training"

# Apply training mask to extract bins for training
bedtools intersect -v -header -wa \
  -a $inbed \
  -b $exbed \
| bgzip -c \
> $prefix.bed.gz
tabix -f $prefix.bed.gz


###########################################
#    Step 3. Train mutation rate model    #
###########################################
# Note: this step requires data pooled across all chromosomes, so filtered BEDs 
# must be downloaded from an existing source (e.g., Cromwell output bucket).
# This step is *not* necessary within the WDL (these files are localized as inputs)
cromwell_output_bucket="gs://dsmap/cromwell_outputs/TrainMuModel/687e362b-3352-42a8-af50-89e244423fd6"
gsutil -m cp ${cromwell_output_bucket}**/$main_prefix.$cnv.pairs.wCounts.*.training.bed.gz ./

# Set parameters as required for TrainModel
export model="logit"
export prefix="$main_prefix.$cnv"
export contigs_fai="hg38.contigs.dev.fai"

# Build list of training BEDs per contig
while read contig; do
  find / -name "*.pairs.wCounts.$contig.training.bed.gz" \
  | awk -v OFS="\t" -v contig="$contig" '{ print contig, $0 }'
done < <( cut -f1 ${contigs_fai} ) \
> training_beds.tsv

# Train model
athena mu-train \
  --training-data training_beds.tsv \
  --model ${model} \
  --model-outfile ${prefix}.${model}.trained.pkl \
  --stats-outfile ${prefix}.${model}.training_stats.tsv \
  --calibration-outfile ${prefix}.${model}.calibration.tsv


############################################################
#    Step 4. Apply mutation rate model to all bin-pairs    #
############################################################
# Using the outputs from steps 1 & 3, generated above
# Note: in WDL, the input pairs file will be sharded for improved parallelization
# In practice, we subset to a single shard here for dev. purposes

# Setp parameters as required for ShardAndPredictMu.wdl
prefix="$main_prefix.$cnv"
bed="$prefix.pairs.wCounts.$contig.bed.gz"
shard_size=10000
trained_model="${prefix}.${model}.trained.pkl"

# Shard BED (note that in practice this will be automated by ShardAndPredictMu.wdl)
zcat ${bed} | head -n $(( ${shard_size} + 1 )) | bgzip -c > shard.bed.gz
tabix -f shard.bed.gz

# Predict mutation rates for all bins in shard
athena mu-predict \
  --trained-model ${trained_model} \
  --outfile ${prefix}.mu.bed.gz \
  --bgzip \
  shard.bed.gz








