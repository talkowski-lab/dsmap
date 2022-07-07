# Processing of individual WGS cohorts  

Each WGS-based SV callset was subject to additional curation, quality control, and formatting steps prior to running the DSMap model.  

Summary metadata for each cohort is provided in [the `data/` subdirectory](https://github.com/talkowski-lab/dsmap/tree/main/data).  

This subdirectory contains the code required to process each WGS cohort. Brief notes for each cohort are provided below.  

### HGSV  

Sample exclusion criteria:  
*  Retained 3,202 samples from five superpopulations (AFR, AMR, EAS, EUR, SAS) with at least 100 samples per population  
*  Excluded 612 children from parent-child trios, retaining 2,590 unrelated samples

Variant filtering criteria:
*  Biallelic deletions or duplications  
*  Autosomal  
*  VCF `FILTER` = `PASS`
*  Maximum 1% allele frequency in each of the five superpopulations  
*  Non-null genotypes present in >90% of all samples (>2,331/2,590)
*  Minimum `QUAL` score of 2  
*  Appears in Hardy-Weinberg Equilibrium (_P_ > 0.01 from H-W chi-squared test)  

