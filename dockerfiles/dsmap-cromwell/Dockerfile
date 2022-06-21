# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile (with Cromwell/Cromshell enabled)
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License


##################################
# BUILD STAGE 1: CONDA ENVIRONMENT
##################################

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge AS build
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Load Athena & DSMap repos (must be present in build context)
ADD athena /opt/athena
ADD dsmap /opt/dsmap
RUN cp /opt/athena/config/condarc ~/.condarc

# # Create default python environment
RUN mamba update --all && \
    mamba env create --file=/opt/dsmap/config/environment.dsmap-cromwell.yml

# Use conda-pack to create a standalone enviornment
RUN mamba install conda-pack && \
    conda-pack -n dsmap-cromwell -o /tmp/env.tar && \
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
COPY --from=build /opt/dsmap /opt/dsmap

# Make bash default shell & load DSMap environment on container startup
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN echo "source /venv/bin/activate" >> ~/.bashrc
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install make, zlib, and other utils required for Cromwell/womtool
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      less \
      default-libmysqlclient-dev \
      zip unzip \
      mailutils \
      bsdmainutils \
      curl

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
ARG DEBIAN_FRONTEND=noninteractive
RUN /opt/athena/config/install_gcloud.sh

# Set default Cromwell server URL
ENV CROMWELL_URL https://cromwell-talkowski.dsde-methods.broadinstitute.org/
RUN echo -e "$CROMWELL_URL\nYes" | cromshell

# Cromshell conflicts with Anaconda's version of curl
# Solution: symlink /usr/bin/curl to /opt/bin/curl and add /opt/bin to $PATH
# /opt/bin must be added to path at the end of ~/.bashrc (rather than via ENV PATH)
# due to conda activate prepending to the path when the container is deployed
RUN mkdir /opt/bin && \
    ln -s /usr/bin/curl /opt/bin/curl && \
    echo 'export PATH="/opt/bin:$PATH"' >> ~/.bashrc

# Cromshell within docker is unable to get local issuer certificate
# Can work around this issue by replacing conda-installed cromshell with modified 
# cromshell script that aliases curl -k and allows insecure connections
RUN ln -sf /opt/dsmap/cromwell/cromshell_insecure_curl/cromshell $( which cromshell )

# Launch bash
CMD ["/bin/bash"]
