#!/usr/bin/env bash

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>

# Code to preprocess reference data


################
### CENTROMERES
################
# Download from UCSC browser track
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz

# Create BED of intervals with smallest and largest centromere coordinates per chromosome
zcat < centromeres.txt.gz \
| awk -F'\t' -v OFS='\t' '{print $2, $3, $4}' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools groupby -i hg38.centromeres.all.bed.gz -g 1 -c 2,2 -o min,max \
| bgzip -c \
> hg38.centromeres.bed.gz

# Copy to GCP bucket
gsutil cp hg38.centromeres.bed.gz gs://dsmap/data/references/