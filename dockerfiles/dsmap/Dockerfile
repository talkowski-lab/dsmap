# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: us.gcr.io/broad-dsmap/athena-cloud:latest
ARG ATHENA_CLOUD_BASE_IMAGE="us.gcr.io/broad-dsmap/athena-cloud:latest"
FROM ${ATHENA_CLOUD_BASE_IMAGE}
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install additional libraries
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      zip unzip

# Add DSMap code (assumes repo has been cloned into build context)
ADD dsmap /opt/dsmap

# Launch bash
CMD ["/bin/bash"]
