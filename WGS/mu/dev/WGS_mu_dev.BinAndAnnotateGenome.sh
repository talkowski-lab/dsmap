#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Development code for prototype WGS mutation rate model
# Genome binning and annotation


###############
#    Setup    #
###############
# Launch Docker
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-cromwell:wgs-mu-dev-7411a1-722506

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


########################################
#    Get distributions of CNV data     #
########################################
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
  "BinAndAnnotateGenome.bins_per_pair_shard": 140,
  "BinAndAnnotateGenome.bins_per_shard": 5000,
  "BinAndAnnotateGenome.bin_annotations_list_remote": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_remote.tsv",
  "BinAndAnnotateGenome.bin_annotations_list_ucsc": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_ucsc.tsv",
  "BinAndAnnotateGenome.bin_exclusion_mask": "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.bin_size": 20000,
  "BinAndAnnotateGenome.contigs_fai": "gs://dsmap/data/dev/hg38.contigs.dev.fai",
  "BinAndAnnotateGenome.decompose_features": true,
  "BinAndAnnotateGenome.feature_transformations_tsv": "gs://dsmap/data/dev/WGS.mu.dev.feature_transformations.tsv",
  "BinAndAnnotateGenome.max_pair_distance": 600000,
  "BinAndAnnotateGenome.max_pcs": 1000,
  "BinAndAnnotateGenome.pairs_for_pca": 100000,
  "BinAndAnnotateGenome.pairs_per_shard_apply_pca": 50000,
  "BinAndAnnotateGenome.pair_annotations_list_remote": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_remote.pairs.tsv",
  "BinAndAnnotateGenome.pair_annotations_list_ucsc": "gs://dsmap/data/dev/WGS.mu.dev.athena_tracklist_ucsc.pairs.tsv",
  "BinAndAnnotateGenome.pair_exclusion_mask": "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.pca_min_variance": 0.999,
  "BinAndAnnotateGenome.prefix": "$prefix",
  "BinAndAnnotateGenome.ref_build": "hg38",
  "BinAndAnnotateGenome.ref_fasta": "gs://dsmap/data/references/hg38.fa",
  "BinAndAnnotateGenome.snv_mutrates_tsv": "gs://dsmap/data/resources/snv_mutation_rates.Samocha_2014.tsv.gz",
  "BinAndAnnotateGenome.visualize_features_before_pca": true
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

# Move all final outputs to dsmap gs:// bucket for permanent storage
cromshell -t 20 metadata \
  $( cat /cromwell/logs/$prefix.BinAndAnnotateGenome.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .outputs | sed 's/\"/\n/g' | fgrep "gs://" \
| gsutil -m cp -I gs://dsmap/data/dev/$prefix.BinAndAnnotateGenome/


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
  gs://dsmap/data/dev/WGS.mu.dev.feature_transformations.tsv \
  ./

# Set parameters as required in BinAndAnnotateGenome.wdl
export contigs_fai=hg38.contigs.dev.fai
export bin_exclusion_mask=hg38.nmask.bed.gz
export pair_exclusion_mask=hg38.nmask.bed.gz
export bin_size=20000
export max_pair_distance=600000
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
export feature_transformations_tsv=WGS.mu.dev.feature_transformations.tsv
export pca_min_variance=0.999
export pairs_for_pca=10000


####################################
#    Step 1. Create all 1D bins    #
####################################
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


#####################################################################################
#    Prior to Step 2, determine number of pairs to sample per chromosome for PCA    #
#####################################################################################
total_bp=$( awk -v FS="\t" '{ sum+=$2 }END{ print sum }' ${contigs_fai} )
bp_per_pair=$(( ( $total_bp + ${pairs_for_pca} - 1 ) / ${pairs_for_pca} ))
awk -v bp_per_pair="$bp_per_pair" -v FS="\t" \
  '{ printf "%s\t%s\t%0.0f\n", $1, $2, ($2 + bp_per_pair - 1) / bp_per_pair }' \
  ${contigs_fai} \
> ${prefix}.updated.fai


#############################################
#    Step 2. annotate 1D bins per contig    #
#############################################
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


#####################################################
#    Step 3. Pair 2D bins and add 2D annotations    #
#####################################################
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
  cat ${pair_annotations_list} > local_tracks.pairs.tsv
fi
if [ ! -z ${pair_annotations_list_remote} ]; then
  if [ ! -z ${ref_fasta} ]; then
    remote_options="--ref-fasta ${ref_fasta}"
  else
    remote_options=""
  fi
  athena slice-remote $remote_options \
    --updated-tsv local_slices.pairs.tsv \
    ${pair_annotations_list_remote} \
    regions.bed
  cat local_slices.pairs.tsv >> local_tracks.pairs.tsv
fi

# Build options for athena annotate-pairs
# Note: given that homology searching is computationally intensive, for local
# testing purposes we need to downsample to ~100 random pairs
tabix -H ${prefix}.pairs.bed.gz > ${prefix}.pairs.subset.bed
zcat ${prefix}.pairs.bed.gz | fgrep -v "#" \
| shuf | head -n100 \
| sort -Vrk1,1 -k2,2n -k3,3n \
>> ${prefix}.pairs.subset.bed
bgzip -f ${prefix}.pairs.subset.bed
tabix -f ${prefix}.pairs.subset.bed.gz
export pairs="${prefix}.pairs.subset.bed.gz"
athena_options=""
if [ ! -z local_tracks.pairs.tsv ]; then
  athena_options="$athena_options --track-list local_tracks.pairs.tsv"
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
eval $athena_cmd 2> >(fgrep -v "is marked as paired, but its mate does not occur" >&2)
tabix -f ${prefix}.annotated.pairs.bed.gz


###########################################
#    Step 4. Decompose all annotations    #
###########################################
# [Optional] Visualize distributions of all features prior to PCA
mkdir ${prefix}_feature_hists
athena feature-hists \
  --ignore-columns 3 \
  --transformations-tsv ${feature_transformations_tsv} \
  ${prefix}.annotated.pairs.bed.gz \
  ${prefix}_feature_hists/${prefix}

# Learn PCA transformation from subsampled data
# Note: for local dev purposes, we are just using the output from annotate-pairs, above
# In practice, we would semi-randomly sample pairs from every chromosome and learn
# PCA transformation on those bins
athena eigen-bins \
  -o ${prefix}.annotated.pairs.eigen_train.bed.gz \
  --parameters-outfile ${prefix}.pca_model.pickle \
  --min-variance ${pca_min_variance} \
  --transformations-tsv ${feature_transformations_tsv} \
  --whiten \
  --fill-missing mean \
  --max-components 1000 \
  --stats ${prefix}.pca_stats.txt \
  --bgzip \
  ${prefix}.annotated.pairs.bed.gz


########################################
#    Step 5. Apply PCA to all pairs    #
########################################
# Note: for local dev purposes, we are reapplying the same model from step 4
# to the same pairs from step 4. In practice, this model would be applied to all
# bins, including those not used in PCA fitting
athena eigen-bins \
  -o ${prefix}.annotated.pairs.eigen.bed.gz \
  --precomputed-parameters ${prefix}.pca_model.pickle \
  --bgzip \
  ${prefix}.annotated.pairs.bed.gz

# For local dev purposes only, can confirm that the trained model .pickle was applied
# correctly by checking for any differences between -o outputs from Steps 4 & 5
# (Note that in practice we do not need to generate the -o output at Step 4)
diff ${prefix}.annotated.pairs.eigen_train.bed.gz ${prefix}.annotated.pairs.eigen.bed.gz

# For local dev purposes, we can visualize the decomposed features to assess their scaling
if ! [ -e ${prefix}_feature_hists ]; then
  mkdir ${prefix}_feature_hists
fi
athena feature-hists \
  --ignore-columns 3 \
  ${prefix}.annotated.pairs.eigen.bed.gz \
  ${prefix}_feature_hists/${prefix}.decomped




