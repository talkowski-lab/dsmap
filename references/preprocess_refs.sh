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

curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz

zcat centromeres.txt.gz \
| awk -F'\t' -v OFS='\t' '{print $2, $3, $4}' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> hg38.centromeres.all.bed.gz

# Create BED of intervals with smallest and largest centromere coordinates per chromosome
paste \
    <(bedtools groupby -i hg38.centromeres.all.bed.gz -g 1 -c 2 -o min | cut -f1,2) \
    <(bedtools groupby -i hg38.centromeres.all.bed.gz -g 1 -c 2 -o max | cut -f2) \
| bgzip -c \
> hg38.centromeres.bed.gz

gsutil cp hg38.centromeres.bed.gz gs://dsmap/data/references/