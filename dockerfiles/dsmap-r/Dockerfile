# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile (with R installed and R libraries needed fpr plotting])
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: us.gcr.io/broad-dsmap/dsmap
ARG DSMAP_BASE_IMAGE="us.gcr.io/broad-dsmap/dsmap:latest"
FROM ${DSMAP_BASE_IMAGE}
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Create alternative python environment specifically for R
# This is necessary due to the massive number of new packages in r-essentials
# causing long "solving environment" delays with conda install
RUN conda env create --file=/opt/dsmap/conda/dsmap-r.environment.yml
ENV CONDA_DEFAULT_ENV dsmap-r
RUN sed --in-place 's/conda activate dsmap/conda activate dsmap-r/g' ~/.bashrc
ENV PATH /opt/conda/envs/dsmap-r/bin:$PATH
RUN conda activate dsmap-r

# Install dsmapR package from source
ARG DSMAPR_VERSION=0.1.0
RUN Rscript -e "install.packages('/opt/dsmap/dsmapR/dsmapR_${DSMAPR_VERSION}.tar.gz', type='source', repos=NULL)"

# Install system libraries for plotting
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install libxrender1 libcairo2-dev libxt-dev libx11-dev 
RUN apt-get -qqy upgrade

# Install Microsoft TrueType fonts (including Arial). Necessary for rendering .pngs
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install software-properties-common && \
    add-apt-repository contrib
RUN apt-get -qqy update && \
    apt-get -qqy install ttf-mscorefonts-installer

# Launch bash
CMD ["/bin/bash"]
