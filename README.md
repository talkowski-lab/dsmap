# DSMap
### Genome-wide maps of human copy number mutation rates and dosage sensitivity

Copyright (c) 2021-Present, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  

#### _Note: this repository is under active development_

---  

## Synopsis    

This repository contains the working code and scripts used to:  
1. Aggregate copy-number variant (CNV) data across technologies;
2. Build genome-wide CNV mutation rate maps; and
3. Compute genome-wide dosage sensitivity maps  

---  

## Table of Contents  

| Directory | Description |  
| :--- | :--- |  
| [`annotations/`](https://github.com/talkowski-lab/dsmap/tree/main/annotations) | Processing of genome annotations |  
| [`cma/`](https://github.com/talkowski-lab/dsmap/tree/main/cma) | Processing & analysis of chromosomal microarray datasets |  
| [`cromwell/`](https://github.com/talkowski-lab/dsmap/tree/main/cromwell) | Resources for executing cloud workflows on Cromwell |  
| [`data/`](https://github.com/talkowski-lab/dsmap/tree/main/data) | Manifests of datasets used in this project |  
| [`dockerfiles/`](https://github.com/talkowski-lab/dsmap/tree/main/docker) | Build files for Docker Images used in this project |  
| [`integration/`](https://github.com/talkowski-lab/dsmap/tree/main/integration) | Cross-technology integration of dosage sensitivity maps |  
| [`scripts/`](https://github.com/talkowski-lab/dsmap/tree/main/scripts) | Stand-alone helper scripts |  
| [`wdls/`](https://github.com/talkowski-lab/dsmap/tree/main/wdls) | Stand-alone WDL workflows |  
| [`wes/`](https://github.com/talkowski-lab/dsmap/tree/main/wes) | Processing & analysis of whole-exome sequencing datasets |  
| [`wgs/`](https://github.com/talkowski-lab/dsmap/tree/main/wgs) | Processing & analysis of whole-genome sequencing datasets |  

---  

## Key Dependencies  

#### Athena

This repository relies heavily on the functions provided by [`athena`](https://github.com/talkowski-lab/athena), a command-line toolkit designed to support the DSMap project.  

Please [refer to the `athena` repository](https://github.com/talkowski-lab/athena) for more information on `athena`, its functions, or the Docker Image provided for its deployment. 

#### Docker  

Every analysis in this repository is intended to be run in a Docker container.  

If you do not have Docker installed on your machine, it is recommended to do so prior to running any code contained herein.  

#### WDL/Cromwell  

In many cases, analyses were parallelized in Google Cloud using the [Workflow Description Language (WDL)](https://openwdl.org/) and the [Cromwell Workflow Management System](https://cromwell.readthedocs.io/en/stable/).  

---  
