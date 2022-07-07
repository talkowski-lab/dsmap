# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License


##################################
# BUILD STAGE 1: CONDA ENVIRONMENT
##################################

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge:4.12.0-2 AS build
LABEL stage=build
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Load athena repo (must be present in build context)
ADD athena /opt/athena
COPY athena/config/condarc ~/.condarc

# Create default python environment
RUN mamba update --all && \
    mamba env create --file=/opt/athena/config/environment.yml

# Use conda-pack to create a standalone enviornment
RUN mamba install conda-pack && \
    conda-pack -n olympus -o /tmp/env.tar && \
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

# Make bash default shell & load Olympus on container startup
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN echo "source /venv/bin/activate" >> ~/.bashrc
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install make and zlib
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
<<<<<<< HEAD
      zlib1g-dev \
      less \
      default-libmysqlclient-dev \
      curl \
      wget
=======
      zlib1g-dev
      # libcurl4-openssl-dev

# Create python environment
# Note: in practice, it's more sane to install packages one at a time
# due to intermittent HTTPS timeouts
RUN conda update -n base conda
RUN conda create --name=dsmap python=3.7.3
# RUN conda config --set ssl_verify no
RUN conda config --add channels anaconda
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
ENV CONDA_DEFAULT_ENV dsmap
RUN echo "source activate dsmap" >> ~/.bashrc
ENV PATH /opt/conda/envs/dsmap/bin:$PATH
RUN conda install -n dsmap -y openssl curl
RUN conda install -n dsmap -y numpy scipy pandas scikit-learn
RUN conda install -n dsmap -y matplotlib
RUN conda install -n dsmap -y htslib=1.11 samtools=1.11 pysam=0.16
RUN conda install -n dsmap -y pybedtools=0.8.1 bedtools=2.30.0
RUN conda install -n dsmap -y rapidfuzz
RUN conda install -n dsmap -y bzip2
RUN conda install -n dsmap -y click
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing default-libmysqlclient-dev
RUN conda install -n dsmap -y mysqlclient
RUN conda install -n dsmap -y pybigwig=0.3.17
>>>>>>> main

# Increase max MySQL packet size for large UCSC queries
RUN echo "[mysqld]\nmax_allowed_packet = 2G" > /etc/my.cnf && \
    echo "[client]\nmax_allowed_packet = 2G" >> /etc/my.cnf

# Developer aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Install Athena
RUN cd /opt/athena && \
    pip install -e .

# Launch bash
CMD ["/bin/bash"]
