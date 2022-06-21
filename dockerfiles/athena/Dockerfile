# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License


##################################
# BUILD STAGE 1: CONDA ENVIRONMENT
##################################

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge AS build
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
FROM debian:bullseye AS runtime
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
      zlib1g-dev \
      less \
      default-libmysqlclient-dev

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
