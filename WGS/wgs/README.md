# Whole-genome sequencing data manifest  

All WGS-related data is stored in the following Google Bucket:  
`gs://dsmap/data/WGS/`  

Additional details for each cohort are provided in the table below:  

| Cohort | Abbreviation | Used in Final Model | Samples | Description | Citation | Original Source | Notes |    
| :--- | :--- | :--- | ---: | :--- | :--- | :--- | :--- |  
| Human Genome Structural Variation Consortium | HGSV | No | TBD | GATK-SV callset on N,NNN samples from the extended 1000 Genomes Project resequenced to deep (~30X) coverage at the New York Genome Center | [Byrska-Bishop _et al_., _bioRxiv_ (2021)](https://doi.org/10.1101/2021.02.06.430068) | `gs://talkowski-sv-gnomad-output/1KGP/final_vcf/1KGP_2504_and_698_with_GIAB.boost.PASS_gt_revised.vcf.gz` | Dataset used for pipeline development only; not included in final DSMap model |  
