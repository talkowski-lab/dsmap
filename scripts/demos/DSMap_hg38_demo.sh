#!/usr/bin/env bash

#########################
#    DSMap hg38 Demo    #
#########################

# Copyright (c) 2023-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rcollins@broadinstitute.org>

# DSMap hg38 demonstration project
# gnomAD-SV v3.0 


# Set local parameters
export prefix=DSMap_hg38_demo
export WRKDIR=/Users/collins/Desktop/Collins/Talkowski/NGS/SV_Projects/dosage_sensitivity/$prefix
cd $WRKDIR
export DSMAP=/Users/collins/Desktop/Collins/Talkowski/NGS/SV_Projects/dosage_sensitivity/dsmap
export GATKSV=/Users/collins/Desktop/Collins/Talkowski/code/gatk-sv


# Prep local directory structure
for subdir in cromwell cromwell/logs cromwell/outputs cromwell/inputs scratch; do
  mkdir ${WRKDIR}/$subdir
done


#################################
### STEP 1: FILTER GNOMAD VCF ###
#################################
# Make inputs .json for creating sites VCFs
cat <<EOF > ${WRKDIR}/cromwell/inputs/MakeSitesVcf.input.json
{
  "MakeSitesVcf.athena_docker":
    "us.gcr.io/broad-dsmap/athena:wgs-ds-dev-69eab0-48154f",
  "MakeSitesVcf.vcfs": 
    [
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-0/cacheCopy/gnomAD_SV_v3.releasable.chr1.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-1/cacheCopy/gnomAD_SV_v3.releasable.chr2.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-2/cacheCopy/gnomAD_SV_v3.releasable.chr3.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-3/cacheCopy/gnomAD_SV_v3.releasable.chr4.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-4/cacheCopy/gnomAD_SV_v3.releasable.chr5.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-5/cacheCopy/gnomAD_SV_v3.releasable.chr6.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-6/cacheCopy/gnomAD_SV_v3.releasable.chr7.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-7/cacheCopy/gnomAD_SV_v3.releasable.chr8.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-8/cacheCopy/gnomAD_SV_v3.releasable.chr9.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-9/cacheCopy/gnomAD_SV_v3.releasable.chr10.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-10/cacheCopy/gnomAD_SV_v3.releasable.chr11.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-11/cacheCopy/gnomAD_SV_v3.releasable.chr12.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-12/cacheCopy/gnomAD_SV_v3.releasable.chr13.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-13/cacheCopy/gnomAD_SV_v3.releasable.chr14.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-14/cacheCopy/gnomAD_SV_v3.releasable.chr15.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-15/cacheCopy/gnomAD_SV_v3.releasable.chr16.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-16/cacheCopy/gnomAD_SV_v3.releasable.chr17.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-17/cacheCopy/gnomAD_SV_v3.releasable.chr18.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-18/cacheCopy/gnomAD_SV_v3.releasable.chr19.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-19/cacheCopy/gnomAD_SV_v3.releasable.chr20.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-20/cacheCopy/gnomAD_SV_v3.releasable.chr21.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-21/cacheCopy/gnomAD_SV_v3.releasable.chr22.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-22/cacheCopy/gnomAD_SV_v3.releasable.chrX.annotated.header_fix.MEI_DEL_patch.vcf.gz",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-23/cacheCopy/gnomAD_SV_v3.releasable.chrY.annotated.header_fix.MEI_DEL_patch.vcf.gz"
    ],
  "MakeSitesVcf.vcf_idxs": 
    [
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-0/cacheCopy/gnomAD_SV_v3.releasable.chr1.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-1/cacheCopy/gnomAD_SV_v3.releasable.chr2.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-2/cacheCopy/gnomAD_SV_v3.releasable.chr3.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-3/cacheCopy/gnomAD_SV_v3.releasable.chr4.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-4/cacheCopy/gnomAD_SV_v3.releasable.chr5.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-5/cacheCopy/gnomAD_SV_v3.releasable.chr6.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-6/cacheCopy/gnomAD_SV_v3.releasable.chr7.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-7/cacheCopy/gnomAD_SV_v3.releasable.chr8.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-8/cacheCopy/gnomAD_SV_v3.releasable.chr9.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-9/cacheCopy/gnomAD_SV_v3.releasable.chr10.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-10/cacheCopy/gnomAD_SV_v3.releasable.chr11.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-11/cacheCopy/gnomAD_SV_v3.releasable.chr12.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-12/cacheCopy/gnomAD_SV_v3.releasable.chr13.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-13/cacheCopy/gnomAD_SV_v3.releasable.chr14.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-14/cacheCopy/gnomAD_SV_v3.releasable.chr15.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-15/cacheCopy/gnomAD_SV_v3.releasable.chr16.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-16/cacheCopy/gnomAD_SV_v3.releasable.chr17.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-17/cacheCopy/gnomAD_SV_v3.releasable.chr18.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-18/cacheCopy/gnomAD_SV_v3.releasable.chr19.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-19/cacheCopy/gnomAD_SV_v3.releasable.chr20.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-20/cacheCopy/gnomAD_SV_v3.releasable.chr21.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-21/cacheCopy/gnomAD_SV_v3.releasable.chr22.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-22/cacheCopy/gnomAD_SV_v3.releasable.chrX.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi",
      "gs://talkowski-sv-gnomad-output/zero/AnnotateVcf/07ebdac2-d3c1-451f-8b26-3d322ad462a0/call-MeiDelPatch/shard-23/cacheCopy/gnomAD_SV_v3.releasable.chrY.annotated.header_fix.MEI_DEL_patch.vcf.gz.tbi"
    ]
}
EOF

# Make dependencies .zip
if [ -e ${WRKDIR}/dsmap.dependencies.zip ]; then
  rm ${WRKDIR}/dsmap.dependencies.zip
fi
cd ${DSMAP}/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip ${WRKDIR}/ && \
cd -

# Submit to cromwell
cromshell -t 20 submit \
  ${DSMAP}/wdls/MakeSitesVcf.wdl \
  ${WRKDIR}/cromwell/inputs/MakeSitesVcf.input.json \
  ${DSMAP}/cromwell/dsmap.cromwell_options.json \
  ${WRKDIR}/dsmap.dependencies.zip \
> ${WRKDIR}/cromwell/logs/MakeSitesVcf.log

# Check status
cromshell -t 20 metadata \
  $( cat ${WRKDIR}/cromwell/logs/MakeSitesVcf.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status

# Download all sites.vcfs, merge locally, and upload merged VCF to google bucket
mkdir $WRKDIR/scratch/site_vcfs_premerge/
cromshell -t 20 metadata \
  $( cat ${WRKDIR}/cromwell/logs/MakeSitesVcf.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .outputs | fgrep "gs://" | gtr -d '",' | gsed 's/\ //g' \
| gsutil -m cp -I $WRKDIR/scratch/site_vcfs_premerge/
bcftools concat --naive \
  -O z -o $WRKDIR/scratch/gnomad_v3.sv.sites.vcf.gz \
  $WRKDIR/scratch/site_vcfs_premerge/*vcf.gz
tabix -p vcf -f $WRKDIR/scratch/gnomad_v3.sv.sites.vcf.gz
gsutil -m cp \
  $WRKDIR/scratch/gnomad_v3.sv.sites.vcf.gz* \
  gs://dsmap/hg38_demo/data/WGS/
rm -rf \
  $WRKDIR/scratch/site_vcfs_premerge \
  $WRKDIR/scratch/gnomad_v3.sv.sites.vcf.gz*

# Setup for filtering
docker run --rm -it us.gcr.io/broad-dsmap/athena-cloud:wgs-mu-dev-hg38-97e064-ae06a9
gcloud auth login

# Download sites VCF
gsutil -m cp \
  gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.sites.vcf.gz* \
  ./

# Filter
for CNV in DEL DUP; do
  athena vcf-filter \
    --include-chroms "$( seq 1 22 | awk '{ print "chr"$0 }' | paste -s -d, )" \
    --svtypes $CNV \
    --maxAF 0.01 \
    --af-field afr_AF \
    --af-field amr_AF \
    --af-field eas_AF \
    --af-field nfe_AF \
    --af-field sas_AF \
    --minAC 1 \
    --minAN 56690 \
    --minQUAL 2 \
    --pHWE 0.000001 \
    --bgzip \
    gnomad_v3.sv.sites.vcf.gz \
    gnomad_v3.sv.filtered.$CNV.sites.vcf.gz
  tabix -f gnomad_v3.sv.filtered.$CNV.sites.vcf.gz
done

# Copy to GCP bucket
gsutil -m cp \
  gnomad_v3.sv.filtered.*.sites.vcf.gz* \
  gs://dsmap/hg38_demo/data/WGS/

# Get stats
for CNV in DEL DUP; do
  echo -e "\n\n$CNV"
  athena vcf-stats gnomad_v3.sv.filtered.$CNV.sites.vcf.gz
done


######################################
### STEP 2: PREPROCESS ENCODE DATA ###
######################################
# Write list of ENCODE tracks to process
for wrapper in 1; do
  echo -e "https://www.encodeproject.org/files/ENCFF254YHN/@@download/ENCFF254YHN.bam\tENCFF254YHN_testis_total_RNAseq"
  echo -e "https://www.encodeproject.org/files/ENCFF177EPU/@@download/ENCFF177EPU.bam\tENCFF177EPU_ovary_total_RNAseq"
  echo -e "https://www.encodeproject.org/files/ENCFF423JSR/@@download/ENCFF423JSR.bam\tENCFF423JSR_testis_ATACseq"
  echo -e "https://www.encodeproject.org/files/ENCFF053RTV/@@download/ENCFF053RTV.bam\tENCFF053RTV_ovary_ATACseq"
done | awk -v OFS="\t" '{ print $0, "gs://dsmap/hg38_demo/data/ENCODE/bams/" }' \
> ${WRKDIR}/scratch/ENCODE.bams_to_curate.tsv
gsutil cp \
  ${WRKDIR}/scratch/ENCODE.bams_to_curate.tsv \
  gs://dsmap/hg38_demo/data/ENCODE/bams/

# Make inputs .json
cat <<EOF > ${WRKDIR}/cromwell/inputs/ENCODE.PreprocessExternalBAMs.input.json
{
  "PreprocessExternalBAMs.athena_cloud_docker": "us.gcr.io/broad-dsmap/athena-cloud:latest",
  "PreprocessExternalBAMs.bam_urls_tsv": "gs://dsmap/hg38_demo/data/ENCODE/bams/ENCODE.bams_to_curate.tsv"
}
EOF

# Make dependencies .zip
if [ -e ${WRKDIR}/dsmap.dependencies.zip ]; then
  rm ${WRKDIR}/dsmap.dependencies.zip
fi
cd ${DSMAP}/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip ${WRKDIR}/ && \
cd -

# Submit to cromwell
cromshell -t 20 submit \
  ${DSMAP}/wdls/PreprocessExternalBAMs.wdl \
  ${WRKDIR}/cromwell/inputs/ENCODE.PreprocessExternalBAMs.input.json \
  ${DSMAP}/cromwell/dsmap.cromwell_options.json \
  ${WRKDIR}/dsmap.dependencies.zip \
> ${WRKDIR}/cromwell/logs/ENCODE.PreprocessExternalBAMs.log

# Check status
cromshell -t 20 metadata \
  $( cat ${WRKDIR}/cromwell/logs/ENCODE.PreprocessExternalBAMs.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status


#######################################
### STEP 3: BIN AND ANNOTATE GENOME ###
#######################################
# Make inputs .json
cat <<EOF > ${WRKDIR}/cromwell/inputs/$prefix.BinAndAnnotateGenome.input.json
{
  "BinAndAnnotateGenome.athena_cloud_docker": 
    "us.gcr.io/broad-dsmap/athena-cloud:wgs-ds-dev-69eab0-48154f",
  "BinAndAnnotateGenome.athena_docker": 
    "us.gcr.io/broad-dsmap/athena:wgs-ds-dev-f0740e-48154f",
  "BinAndAnnotateGenome.bins_per_pair_shard": 
    150,
  "BinAndAnnotateGenome.bins_per_shard": 
    2500,
  "BinAndAnnotateGenome.bin_annotations_list_remote": 
    "gs://dsmap/hg38_demo/data/tracklists/$prefix.athena_tracklist_remote.tsv",
  "BinAndAnnotateGenome.bin_annotations_list_ucsc": 
    "gs://dsmap/hg38_demo/data/tracklists/$prefix.athena_tracklist_ucsc.tsv",
  "BinAndAnnotateGenome.bin_exclusion_mask": 
    "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.bin_size": 
    5000,
  "BinAndAnnotateGenome.contigs_fai": 
    "gs://dsmap/data/references/hg38.autosomes.fai",
  "BinAndAnnotateGenome.decompose_features": 
    true,
  "BinAndAnnotateGenome.feature_transformations_tsv": 
    "gs://dsmap/hg38_demo/data/tracklists/$prefix.feature_transformations.tsv",
  "BinAndAnnotateGenome.max_pair_distance": 
    500000,
  "BinAndAnnotateGenome.max_pcs": 
    1000,
  "BinAndAnnotateGenome.pairs_for_pca": 
    10000000,
  "BinAndAnnotateGenome.pairs_per_shard_apply_pca": 
    50000,
  "BinAndAnnotateGenome.pair_annotations_list_ucsc": 
    "gs://dsmap/hg38_demo/data/tracklists/$prefix.athena_tracklist_ucsc.pairs.tsv",
  "BinAndAnnotateGenome.pair_exclusion_mask": 
    "gs://dsmap/data/references/hg38.nmask.bed.gz",
  "BinAndAnnotateGenome.pca_min_variance": 
    0.999,
  "BinAndAnnotateGenome.prefix": 
    "$prefix",
  "BinAndAnnotateGenome.ref_build": 
    "hg38",
  "BinAndAnnotateGenome.ref_fasta": 
    "gs://dsmap/data/references/hg38.fa",
  "BinAndAnnotateGenome.snv_mutrates_tsv": 
    "gs://dsmap/data/resources/snv_mutation_rates.Samocha_2014.tsv.gz",
  "BinAndAnnotateGenome.visualize_features_before_pca": 
    true,
  "BinAndAnnotateGenome.runtime_attr_learn_pca":
    {"mem_gb" : 32, "cpu_cores" : 8, "disk_gb" : 100},
  "BinAndAnnotateGenome.runtime_attr_visualize_features":
    {"mem_gb" : 32, "cpu_cores" : 8, "disk_gb" : 100},
  "BinAndAnnotateGenome.athena_docker_dev": 
    "us.gcr.io/broad-dsmap/athena:wgs-mu-dev-hg38-72c115-41ba90"
}
EOF

# Make dependencies .zip
if [ -e ${WRKDIR}/dsmap.dependencies.zip ]; then
  rm ${WRKDIR}/dsmap.dependencies.zip
fi
cd ${DSMAP}/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip ${WRKDIR}/ && \
cd -

# Submit to cromwell
cromshell -t 20 submit \
  ${DSMAP}/wdls/BinAndAnnotateGenome.wdl \
  ${WRKDIR}/cromwell/inputs/$prefix.BinAndAnnotateGenome.input.json \
  ${DSMAP}/cromwell/dsmap.cromwell_options.json \
  ${WRKDIR}/dsmap.dependencies.zip \
> ${WRKDIR}/cromwell/logs/$prefix.BinAndAnnotateGenome.log

# Check status
cromshell -t 120 metadata \
  $( cat ${WRKDIR}/cromwell/logs/$prefix.BinAndAnnotateGenome.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status

# Move all final outputs to dsmap gs:// bucket for permanent storage
workflow_id=$( cat ${WRKDIR}/cromwell/logs/$prefix.BinAndAnnotateGenome.log | tail -n1 | jq .id | tr -d '"' )
cromshell -t 20 metadata "$workflow_id" \
| jq .outputs | sed 's/\"/\n/g' | fgrep "gs://" \
| gsutil -m cp -I gs://dsmap/hg38_demo/data/$prefix.BinAndAnnotateGenome/
# Also make note of most recent workflow ID (in order to reference intermediate 
# outputs after Docker container is closed)
echo "$workflow_id" > $prefix.BinAndAnnotateGenome.latest_workflow_id.txt
gsutil -m cp \
  $prefix.BinAndAnnotateGenome.latest_workflow_id.txt \
  gs://dsmap/hg38_demo/data/$prefix.BinAndAnnotateGenome/


#########################################
### STEP 4: TRAIN MUTATION RATE MODEL ###
#########################################
# Make dependencies .zip
if [ -e ${WRKDIR}/dsmap.dependencies.zip ]; then
  rm ${WRKDIR}/dsmap.dependencies.zip
fi
cd ${DSMAP}/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip ${WRKDIR}/ && \
cd -

# Make athena mu-train config
cat << EOF > ${WRKDIR}/$prefix.athena_mu_train.config.preMuPrior.json
{
  "model_class" : "logit",
  "cv_eval" : true,
  "max_cv_k" : 22,
  "l2" : 0.1,
  "max_epochs" : 100000,
  "seed" : 2023,
  "gw_mu_prior" : MU_PRIOR_HERE
}
EOF
gsed 's/MU_PRIOR_HERE/0.146/g' \
  ${WRKDIR}/$prefix.athena_mu_train.config.preMuPrior.json \
> ${WRKDIR}/$prefix.athena_mu_train.DEL.config.json
gsed 's/MU_PRIOR_HERE/0.040/g' \
  ${WRKDIR}/$prefix.athena_mu_train.config.preMuPrior.json \
> ${WRKDIR}/$prefix.athena_mu_train.DUP.config.json
gsutil cp \
  ${WRKDIR}/$prefix.athena_mu_train.*.config.json \
  gs://dsmap/hg38_demo/config/

# Must run once for each CNV type
for cnv in DEL DUP; do
  # Make inputs .json
  cat <<EOF > ${WRKDIR}/cromwell/inputs/$prefix.TrainMuModel.$cnv.input.json
{
  "TrainMuModel.athena_cloud_docker": 
    "us.gcr.io/broad-dsmap/athena-cloud:wgs-mu-dev-hg38-72c115-41ba90",
  "TrainMuModel.athena_docker": 
    "us.gcr.io/broad-dsmap/athena:wgs-mu-dev-hg38-72c115-41ba90",
  "TrainMuModel.athena_training_config":
    "gs://dsmap/hg38_demo/config/$prefix.athena_mu_train.$cnv.config.json",
  "TrainMuModel.cnv": 
    "$cnv",
  "TrainMuModel.contigs_fai": 
    "gs://dsmap/data/references/hg38.autosomes.fai",
  "TrainMuModel.dsmap_r_docker": 
    "us.gcr.io/broad-dsmap/dsmap-r:wgs-ds-dev-f0740e-b3d902",
  "TrainMuModel.model": 
    "logit",
  "TrainMuModel.pairs_bed_prefix": 
    "$prefix.pairs.eigen",
  "TrainMuModel.pairs_bucket": 
    "gs://dsmap/hg38_demo/data/$prefix.BinAndAnnotateGenome",
  "TrainMuModel.prefix": 
    "$prefix",
  "TrainMuModel.run_diagnostics": 
    true,
  "TrainMuModel.shard_size_apply_mu": 
    100000,
  "TrainMuModel.training_mask": 
    "gs://dsmap/data/athena/athena.hg38.training_exclusion.bed.gz",
  "TrainMuModel.vcf": 
    "gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.filtered.$cnv.sites.vcf.gz",
  "TrainMuModel.vcf_idx": 
    "gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.filtered.$cnv.sites.vcf.gz.tbi",
  "TrainMuModel.runtime_attr_intersect_svs":
    {"mem_gb": 16, "cpu_cores" : 4},
  "TrainMuModel.runtime_attr_diagnostics":
    {"mem_gb": 16, "cpu_cores" : 4}
}
EOF

  # Submit to cromwell
  cromshell -t 20 submit \
    ${DSMAP}/wdls/TrainMuModel.wdl \
    ${WRKDIR}/cromwell/inputs/$prefix.TrainMuModel.$cnv.input.json \
    ${DSMAP}/cromwell/dsmap.cromwell_options.json \
    ${WRKDIR}/dsmap.dependencies.zip \
  > ${WRKDIR}/cromwell/logs/$prefix.TrainMuModel.$cnv.log
done

# Check status
for cnv in DEL DUP; do
  cromshell -t 20 metadata \
    $( cat ${WRKDIR}/cromwell/logs/$prefix.TrainMuModel.$cnv.log | tail -n1 | jq .id | tr -d '"' ) \
  | jq .status
done

# Get timing diagrams
for cnv in DEL DUP; do
  cromshell -t 20 timing \
    $( cat ${WRKDIR}/cromwell/logs/$prefix.TrainMuModel.$cnv.log | tail -n1 | jq .id | tr -d '"' )
done

# Move all final outputs to dsmap gs:// bucket for permanent storage
for cnv in DEL DUP; do
  workflow_id=$( cat ${WRKDIR}/cromwell/logs/$prefix.TrainMuModel.$cnv.log | tail -n1 | jq .id | tr -d '"' )
  cromshell -t 20 metadata $workflow_id \
  | jq .outputs | sed 's/\"/\n/g' | fgrep "gs://" \
  | gsutil -m cp -I gs://dsmap/hg38_demo/data/$prefix.TrainMuModel/
  # Also make note of most recent workflow ID (in order to reference intermediate 
  # outputs after Docker container is closed)
  echo "$workflow_id" > $prefix.$cnv.TrainMuModel.latest_workflow_id.txt
  gsutil -m cp \
    $prefix.$cnv.TrainMuModel.latest_workflow_id.txt \
    gs://dsmap/hg38_demo/data/$prefix.TrainMuModel/
done


################################################
### STEP 5: COMPUTE GENIC DOSAGE SENSITIVITY ###
################################################
# Make inputs .json
cat <<EOF > ${WRKDIR}/cromwell/inputs/$prefix.CalcGenicDosageSensitivity.input.json
{
  "CalcGenicDosageSensitivity.athena_docker": 
    "us.gcr.io/broad-dsmap/athena:wgs-mu-dev-hg38-72c115-41ba90",
  "CalcGenicDosageSensitivity.contigs_fai": 
    "gs://dsmap/data/references/hg38.autosomes.fai",
  "CalcGenicDosageSensitivity.del_vcf": 
    "gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.filtered.DEL.sites.vcf.gz",
  "CalcGenicDosageSensitivity.del_vcf_idx": 
    "gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.filtered.DEL.sites.vcf.gz.tbi",
  "CalcGenicDosageSensitivity.dsmap_r_docker": 
    "us.gcr.io/broad-dsmap/dsmap-r:wgs-ds-dev-f0740e-b3d902",
  "CalcGenicDosageSensitivity.dup_vcf": 
    "gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.filtered.DUP.sites.vcf.gz",
  "CalcGenicDosageSensitivity.dup_vcf_idx": 
    "gs://dsmap/hg38_demo/data/WGS/gnomad_v3.sv.filtered.DUP.sites.vcf.gz.tbi",
  "CalcGenicDosageSensitivity.gtf": 
    "gs://dsmap/hg38_demo/data/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf.gz",
  "CalcGenicDosageSensitivity.gtf_idx": 
    "gs://dsmap/hg38_demo/data/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf.gz.tbi",
  "CalcGenicDosageSensitivity.mu_bed_prefix": 
    "$prefix",
  "CalcGenicDosageSensitivity.mu_bucket": 
    "gs://dsmap/hg38_demo/data/$prefix.TrainMuModel",
  "CalcGenicDosageSensitivity.prefix": 
    "$prefix",
  "CalcGenicDosageSensitivity.runtime_attr_query_mu": 
    {"mem_gb" : 32, "cpu_cores" : 8, "boot_disk_gb" : 15}
}
EOF

# Make dependencies .zip
if [ -e ${WRKDIR}/dsmap.dependencies.zip ]; then
  rm ${WRKDIR}/dsmap.dependencies.zip
fi
cd ${DSMAP}/wdls/ && \
zip dsmap.dependencies.zip *wdl && \
mv dsmap.dependencies.zip ${WRKDIR}/ && \
cd -

# Submit to cromwell
cromshell -t 20 submit \
  ${DSMAP}/wdls/CalcGenicDosageSensitivity.wdl \
  ${WRKDIR}/cromwell/inputs/$prefix.CalcGenicDosageSensitivity.input.json \
  ${DSMAP}/cromwell/dsmap.cromwell_options.json \
  ${WRKDIR}/dsmap.dependencies.zip \
> ${WRKDIR}/cromwell/logs/$prefix.CalcGenicDosageSensitivity.log

# Check status
cromshell -t 20 metadata \
  $( cat ${WRKDIR}/cromwell/logs/$prefix.CalcGenicDosageSensitivity.log | tail -n1 | jq .id | tr -d '"' ) \
| jq .status

# Move all final outputs to dsmap gs:// bucket for permanent storage
workflow_id=$( cat ${WRKDIR}/cromwell/logs/$prefix.CalcGenicDosageSensitivity.log | tail -n1 | jq .id | tr -d '"' )
cromshell -t 20 metadata $workflow_id \
| jq .outputs | sed 's/\"/\n/g' | fgrep "gs://" \
| gsutil -m cp -I N
# Also make note of most recent workflow ID (in order to reference intermediate 
# outputs after Docker container is closed)
echo "$workflow_id" > $prefix.CalcGenicDosageSensitivity.latest_workflow_id.txt
gsutil -m cp \
  $prefix.CalcGenicDosageSensitivity.latest_workflow_id.txt \
  gs://dsmap/hg38_demo/data/$prefix.CalcGenicDosageSensitivity/