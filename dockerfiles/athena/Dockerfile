# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Load athena repo (must be present in build context)
ADD athena /opt/athena

# Install make and zlib
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing default-libmysqlclient-dev

# Create default python environment
# Note: given the large number of packages installed at this step, this environment
# is installed from a config file provided in /opt/athena/conda/environment.yml
RUN conda update -n base conda
RUN conda env create --file=/opt/athena/conda/environment.yml
# TODO: install directly into base environment to decrease layer size
# TODO: install hic-straw from GitHub (https://github.com/aidenlab/straw) *NOT* pip
ENV CONDA_DEFAULT_ENV dsmap
RUN echo "conda activate dsmap" >> ~/.bashrc
ENV PATH /opt/conda/envs/dsmap/bin:$PATH

# Increase max MySQL packet size for large UCSC queries
RUN echo "[mysqld]\nmax_allowed_packet = 2G" > /etc/my.cnf
RUN echo "[client]\nmax_allowed_packet = 2G" >> /etc/my.cnf

# Install Athena
RUN cd /opt/athena && \
    pip install -e .

# Developer aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Clean conda caches to decrease image size
RUN conda clean -y --all

# Launch bash
CMD ["/bin/bash"]
