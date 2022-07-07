# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile (with Google Cloud enabled)
# Copyright (c) 2021-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: us.gcr.io/broad-dsmap/athena
ARG ATHENA_BASE_IMAGE="us.gcr.io/broad-dsmap/athena:latest"
FROM ${ATHENA_BASE_IMAGE}
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install gcloud
ARG DEBIAN_FRONTEND=noninteractive
RUN /opt/athena/config/install_gcloud.sh

# Launch bash
CMD ["/bin/bash"]
