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
      zlib1g-dev \
      libcurl4-openssl-dev

# Create python environment
# Note: in practice, it's more sane to install packages one at a time
# due to intermittent HTTPS timeouts
RUN conda update -n base conda
RUN conda create --name=dsmap python=3.7
RUN conda config --set ssl_verify no
ENV CONDA_DEFAULT_ENV dsmap
RUN echo "source activate dsmap" >> ~/.bashrc
ENV PATH /opt/conda/envs/dsmap/bin:$PATH
RUN conda install -n dsmap -y -c bioconda numpy scipy pandas scikit-learn
RUN conda install -n dsmap -y -c bioconda matplotlib
RUN conda install -n dsmap -y -c bioconda samtools=1.9 pysam=0.15.3
RUN conda install -n dsmap -y -c bioconda pybedtools bedtools=2.27.0
RUN conda install -n dsmap -y -c bioconda bzip2
RUN conda install -n dsmap -y -c bioconda click
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing default-libmysqlclient-dev
RUN conda install -n dsmap -y -c bioconda mysqlclient
RUN conda install -n dsmap -y -c bioconda pybigwig=0.3.17

# Increase max MySQL packet size for large UCSC queries
RUN echo "[mysqld]\nmax_allowed_packet = 2G" > /etc/my.cnf
RUN echo "[client]\nmax_allowed_packet = 2G" >> /etc/my.cnf

# Clone & install Athena
ARG ATHENA_COMMIT="master"
RUN git clone https://github.com/talkowski-lab/athena /opt/athena && \
    cd /opt/athena && \
    git checkout ${ATHENA_COMMIT} && \
    cd -
RUN cd /opt/athena && \
    pip install -e .

# Symlink Athena core reference data to home directory (for GitHub tutorial)
RUN ln -s /opt/athena/data /data

# Developer aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]
