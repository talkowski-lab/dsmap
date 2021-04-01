# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile (with Cromwell/Cromshell enabled)
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: us.gcr.io/broad-dsmap/dsmap
ARG DSMAP_BASE_IMAGE="us.gcr.io/broad-dsmap/dsmap:latest"
FROM ${DSMAP_BASE_IMAGE}
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install libraries required for Cromshell
RUN apt-get -qqy clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      zip unzip \
      mailutils \
      bsdmainutils \
      curl

# Install cromshell & set default server URL
RUN conda install -n dsmap -y cromshell
ENV CROMWELL_URL https://cromwell-talkowski.dsde-methods.broadinstitute.org/
RUN echo "$CROMWELL_URL\nYes" | cromshell

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
