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
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-045d9d-0f58ab

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap
gcloud auth application-default login

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
  "BinAndAnnotateGenome.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
  "BinAndAnnotateGenome.prefix": "$prefix",
  "BinAndAnnotateGenome.bin_size": 25000,
  "BinAndAnnotateGenome.bins_per_shard": 1000,
  "BinAndAnnotateGenome.athena_docker": "us.gcr.io/broad-dsmap/athena:wgs-mu-dev",
  "BinAndAnnotateGenome.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:wgs-mu-dev"
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
  gs://dsmap/data/references/hg38.fa \
  gs://dsmap/data/resources/snv_mutation_rates.Samocha_2014.tsv.gz \
  gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_remote.tsv \
  gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_ucsc.tsv \
  ./

# Set parameters as required in BinAndAnnotateGenome.wdl
export contigs_fai=hg38.contigs.dev.fai
export bin_exclusion_mask=hg38.nmask.bed.gz
export bin_size=25000
export bins_per_shard=1000
export query_slop=1000
export ref_fasta=hg38.fa
export ref_build=hg38
export snv_mutrates_tsv=snv_mutation_rates.Samocha_2014.tsv.gz
export bin_annotations_list_remote=WGS.mu.dev.athena_tracklist_remote.tsv
export bin_annotations_list_ucsc=WGS.mu.dev.athena_tracklist_ucsc.tsv

# Step 1. Create all 1D bins
cut -f1-2 ${contigs_fai} > contigs.genome
export bedtools_genome_file=contigs.genome
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
export bed="${prefix}.${contig}.shard_000000.bed.gz"
# Prepare inputs
out_prefix=$( echo "${bed}" | sed 's/\.bed\.gz//g' )
zcat ${bed} \
| fgrep -v "#" \
| bedtools slop -i - -g ${bedtools_genome_file} -b ${query_slop} \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> regions.bed
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
if [ ! -z ${ref_fasta} ]; then
  samtools faidx ${ref_fasta}
fi

# Slice large input files hosted remotely with athena slice-remote
if [ ! -z ${bin_annotations_list} ]; then
  cat ${bin_annotations_list} > local_tracks.tsv
fi
if [ ! -z ${bin_annotations_list_remote} ]; then
  if [ ! -z ${ref_fasta} ]; then
    remote_options="--ref-fasta ${ref_fasta}"
  else
    remote_options=""
  fi
  athena slice-remote $remote_options \
    --updated-tsv local_slices.tsv \
    ${bin_annotations_list_remote} \
    regions.bed
  cat local_slices.tsv >> local_tracks.tsv
fi

# Build options for athena annotate-bins
athena_options=""
if [ ! -z local_tracks.tsv ]; then
  athena_options="$athena_options --track-list local_tracks.tsv"
fi
if [ ! -z ${bin_annotations_list_ucsc} ]; then
  athena_options="$athena_options --ucsc-list ${bin_annotations_list_ucsc}"
fi
if [ ! -z ${ref_fasta} ]; then
  athena_options="$athena_options --fasta ${ref_fasta}"
fi
if [ ! -z ${snv_mutrates_tsv} ]; then
  athena_options="$athena_options --snv-mutrate ${snv_mutrates_tsv}"
fi

# Annotate bins with athena
athena_cmd="athena annotate-bins --ucsc-ref ${ref_build} $athena_options"
athena_cmd="$athena_cmd --no-ucsc-chromsplit --bgzip"
athena_cmd="$athena_cmd ${bed} ${out_prefix}.annotated.bed.gz"
echo -e "Now annotating using command:\n$athena_cmd"
eval $athena_cmd
tabix -f ${out_prefix}.annotated.bed.gz


