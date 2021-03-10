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
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-db052a-f527e6

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
  "BinAndAnnotateGenome.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:latest",
  "BinAndAnnotateGenome.athena_docker": "us.gcr.io/broad-dsmap/athena:latest",
  "BinAndAnnotateGenome.bins_per_pair_shard": 2500,
  "BinAndAnnotateGenome.bins_per_shard": 5000,
  "BinAndAnnotateGenome.bin_annotations_list_remote": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_remote.tsv",
  "BinAndAnnotateGenome.bin_annotations_list_ucsc": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_ucsc.tsv",
  "BinAndAnnotateGenome.bin_exclusion_mask": "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.bin_size": 20000,
  "BinAndAnnotateGenome.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
  "BinAndAnnotateGenome.max_pair_distance": 700000,
  "BinAndAnnotateGenome.pair_annotations_list_ucsc": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_ucsc.pairs.tsv",
  "BinAndAnnotateGenome.pair_exclusion_mask": "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.prefix": "$prefix",
  "BinAndAnnotateGenome.ref_build": "hg38",
  "BinAndAnnotateGenome.ref_fasta": "gs://dsmap/data/references/hg38.fa",
  "BinAndAnnotateGenome.snv_mutrates_tsv": "gs://dsmap/data/resources/snv_mutation_rates.Samocha_2014.tsv.gz",
  "BinAndAnnotateGenome.visualize_features": true
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
  gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_remote.pairs.tsv \
  gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_ucsc.pairs.tsv \
  ./

# Set parameters as required in BinAndAnnotateGenome.wdl
export contigs_fai=hg38.contigs.dev.fai
export bin_exclusion_mask=hg38.nmask.bed.gz
export pair_exclusion_mask=hg38.nmask.bed.gz
export bin_size=20000
export max_pair_distance=700000
export bins_per_shard=1000
export bins_per_pair_shard=250
export query_slop=1000
export ref_fasta=hg38.fa
export ref_build=hg38
export snv_mutrates_tsv=snv_mutation_rates.Samocha_2014.tsv.gz
export bin_annotations_list_remote=WGS.mu.dev.athena_tracklist_remote.tsv
export bin_annotations_list_ucsc=WGS.mu.dev.athena_tracklist_ucsc.tsv
export pair_annotations_list_remote=WGS.mu.dev.athena_tracklist_remote.pairs.tsv
export pair_annotations_list_ucsc=WGS.mu.dev.athena_tracklist_ucsc.pairs.tsv


#####################################
#    Step 1. Create all 1D bins     #
#####################################
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


##############################################
#    Step 2. annotate 1D bins per contig     #
##############################################
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


######################################################
#    Step 3. Pair 2D bins and add 2D annotations     #
######################################################
# (Note: in WDL, entire chromosome of annotated bins would be passed to this step.
#  For dev. purposes, we are restricting to the single shard of bins defined above)
export all_bins="${out_prefix}.annotated.bed.gz"
zcat $all_bins | head -n $(( $bins_per_pair_shard + 1 )) | bgzip -c > ${prefix}.query_bins.bed.gz
tabix -f ${prefix}.query_bins.bed.gz
export query_bins="${prefix}.query_bins.bed.gz"

# Build options for athena pair-bins
athena_options=""
if [ ! -z ${pair_exclusion_mask} ]; then
  athena_options="$athena_options --exclusion-list ${pair_exclusion_mask}"
fi

# Pair bins with athena
athena_cmd="athena pair-bins --bgzip --max-dist ${max_pair_distance}"
athena_cmd="$athena_cmd --annotate-distance --sort-features --annotate-absdiff"
athena_cmd="$athena_cmd $athena_options --bin-superset ${all_bins}"
athena_cmd="$athena_cmd ${query_bins} ${prefix}.pairs.bed.gz"
echo -e "Now pairing bins using command:\n$athena_cmd"
eval $athena_cmd
tabix -f ${prefix}.pairs.bed.gz

# Prepare inputs for athena annotate-pairs
zcat ${prefix}.pairs.bed.gz \
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
if [ ! -z ${pair_annotations_list} ]; then
  cat ${pair_annotations_list} > local_tracks.tsv
fi
if [ ! -z ${pair_annotations_list_remote} ]; then
  if [ ! -z ${ref_fasta} ]; then
    remote_options="--ref-fasta ${ref_fasta}"
  else
    remote_options=""
  fi
  athena slice-remote $remote_options \
    --updated-tsv local_slices.tsv \
    ${pair_annotations_list_remote} \
    regions.bed
  cat local_slices.tsv >> local_tracks.tsv
fi

# Build options for athena annotate-pairs
export pairs="${prefix}.pairs.bed.gz"
athena_options=""
if [ ! -z local_tracks.tsv ]; then
  athena_options="$athena_options --track-list local_tracks.tsv"
fi
if [ ! -z ${pair_annotations_list_ucsc} ]; then
  athena_options="$athena_options --ucsc-list ${pair_annotations_list_ucsc}"
fi
if [ ! -z ${ref_build} ]; then
  athena_options="$athena_options --ucsc-ref ${ref_build}"
fi
if [ ! -z ${ref_fasta} ]; then
  athena_options="$athena_options --fasta ${ref_fasta}"
fi
# if [ ! -z ${bin_size} ]; then
#   athena_options="$athena_options --binsize ${bin_size}"
# fi

# Annotate pairs with athena
athena_cmd="athena annotate-pairs --bgzip --no-ucsc-chromsplit $athena_options"
athena_cmd="$athena_cmd ${pairs} ${prefix}.annotated.pairs.bed.gz"
echo -e "Now pairing bins using command:\n$athena_cmd"
eval $athena_cmd
tabix -f ${prefix}.annotated.pairs.bed.gz



