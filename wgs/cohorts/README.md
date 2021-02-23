# Processing of individual WGS cohorts  

Each WGS-based SV callset was subject to additional curation, quality control, and formatting steps prior to running the DSMap model.  

Summary metadata for each cohort is provided in [the `data/` subdirectory](https://github.com/talkowski-lab/dsmap/tree/main/data).  

This subdirectory contains the code required to process each WGS cohort. Brief notes for each cohort are provided below.  

### HGSV  

Sample exclusion criteria:  
*  Retained 3,202 samples from five superpopulations (AFR, AMR, EAS, EUR, SAS) with at least 100 samples per population  
