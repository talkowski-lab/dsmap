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
export docker_tag="wgs-mu-dev-d6c705-7be0bf"
docker run --rm -it us.gcr.io/broad-dsmap/dsmap-r:$docker_tag

# Authenticate GCP credentials
gcloud auth login
gcloud config set project broad-dsmap

# Prep output directory structure
mkdir toy_plots/


####################################################################
#    Plot example heatmaps for breakpoint uncertainty schematic    #
####################################################################
# Run plotting code (requires no input data)
/opt/dsmap/paper/toy_plots/plot_toy_bkpt_attribution_heatmaps.R \
  toy_plots/breakpoint_uncertainty


############################################################
#    Copy all outputs to dsmap gs:// bucket for storage    #
############################################################
gsutil -m cp -r \
  toy_plots \
  gs://dsmap/dsmap_paper_files/plots/
