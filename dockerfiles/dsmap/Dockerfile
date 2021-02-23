# DSMap: genome-wide maps of human copy number mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: us.gcr.io/broad-dsmap/athena-cloud:latest
ARG ATHENA_CLOUD_BASE_IMAGE="us.gcr.io/broad-dsmap/athena-cloud:latest"
FROM ${ATHENA_CLOUD_BASE_IMAGE}
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Clone DSMap repo into /opt
ARG DSMAP_COMMIT="main"
RUN git clone https://github.com/talkowski-lab/dsmap /opt/dsmap && \
    cd /opt/dsmap && \
    git checkout ${DSMAP_COMMIT} && \
    cd -

# Launch bash
CMD ["/bin/bash"]
