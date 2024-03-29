# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile (with R installed and R libraries needed fpr plotting])
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License


##################################
# BUILD STAGE 1: CONDA ENVIRONMENT
##################################

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge:4.12.0-2 AS build
LABEL stage=build
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Load Athena & DSMap repos (must be present in build context)
ADD athena /opt/athena
ADD dsmap /opt/dsmap
RUN cp /opt/athena/config/condarc ~/.condarc

# # Create default python environment
RUN mamba update --all && \
    mamba env create --file=/opt/dsmap/config/environment.dsmap-r.yml

# Use conda-pack to create a standalone enviornment
RUN mamba install conda-pack && \
    conda-pack -n dsmap-r -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar -xf /tmp/env.tar && \
    rm /tmp/env.tar && \
    /venv/bin/conda-unpack


##############################
# BUILD STAGE 2: RUNTIME IMAGE
##############################

# Base image: debian bullseye (for smaller runtime image)
FROM ubuntu:focal-20220531 AS runtime
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Copy /venv from the build stage
COPY --from=build /venv /venv
COPY --from=build /opt/athena /opt/athena
COPY --from=build /opt/dsmap /opt/dsmap

# Make bash default shell & load DSMap environment on container startup
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN echo "source /venv/bin/activate" >> ~/.bashrc
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install make and zlib
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      less \
      default-libmysqlclient-dev \
      wget

# Increase max MySQL packet size for large UCSC queries
RUN echo "[mysqld]\nmax_allowed_packet = 2G" > /etc/my.cnf && \
    echo "[client]\nmax_allowed_packet = 2G" >> /etc/my.cnf

# Aliases (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Install Athena
RUN cd /opt/athena && \
    pip install -e .

# Install gcloud
RUN /opt/athena/config/install_gcloud.sh

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
    apt-get -qqy install ttf-mscorefonts-installer

# Launch bash
CMD ["/bin/bash"]
