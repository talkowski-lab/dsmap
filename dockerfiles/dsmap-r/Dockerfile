# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile (with R installed and R libraries needed fpr plotting])
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: us.gcr.io/broad-dsmap/dsmap
ARG DSMAP_BASE_IMAGE="us.gcr.io/broad-dsmap/dsmap:latest"
FROM ${DSMAP_BASE_IMAGE}
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install Microsoft TrueType fonts (including Arial). Necessary for rendering .pngs
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install software-properties-common && \
    add-apt-repository contrib
RUN apt-get -qqy update && \
    apt-get -qqy install ttf-mscorefonts-installer

# Install R base
RUN conda install -n dsmap -y r-base=3.6.3 r-essentials

# Install R packages (note: intentionally do this in chunks for less painful re-builds)
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install libxrender1 libcairo2-dev libxt-dev libx11-dev 
RUN apt-get -qqy upgrade
RUN conda install -n dsmap -y r-optparse r-MASS r-KernSmooth
RUN conda install -n dsmap -y r-caTools r-flux
RUN conda install -n dsmap -y r-Hmisc
RUN conda install -n dsmap -y r-funr r-plotrix
RUN conda install -n dsmap -y r-shape r-gtools
RUN conda install -n dsmap -y r-beeswarm r-vioplot

# Install dsmapR package from source
ARG DSMAPR_VERSION=0.1.0
RUN rscript -e "install.packages('/opt/dsmap/dsmapR/dsmapR_${DSMAPR_VERSION}.tar.gz', type='source', repos=NULL)"

# Launch bash
CMD ["/bin/bash"]
