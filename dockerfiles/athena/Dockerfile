# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install make and zlib
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev
      # libcurl4-openssl-dev

# Create python environment
# Note: in practice, it's more sane to install packages one at a time
# due to intermittent HTTPS timeouts
RUN conda update -n base conda
RUN conda create --name=dsmap python=3.8
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
RUN conda install -n dsmap -y bzip2
RUN conda install -n dsmap -y click
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing default-libmysqlclient-dev
RUN conda install -n dsmap -y mysqlclient
RUN conda install -n dsmap -y pybigwig=0.3.17

# Increase max MySQL packet size for large UCSC queries
RUN echo "[mysqld]\nmax_allowed_packet = 2G" > /etc/my.cnf
RUN echo "[client]\nmax_allowed_packet = 2G" >> /etc/my.cnf

# Install Athena (assumes repo has been cloned into build context)
ADD athena /opt/athena
RUN cd /opt/athena && \
    pip install -e .

# Developer aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]
