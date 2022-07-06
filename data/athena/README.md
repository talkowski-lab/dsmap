# Athena toolkit reference data manifest  

All data files curated for various [`athena`](https://github.com/talkowski-lab/athena) commands throughout these analyses are stored in the following Google Bucket:  
`gs://dsmap/data/athena/`  

Additional details for each dataset are provided in the table below:  

| Dataset | Filename | Citation | Description |  
| :--- | :--- | :--- | :--- |  
| Athena GRCh37 neutral region mask | `athena.GRCh37.training_exclusion.bed.gz` | _N/A_ | A list of regions under purifying selection to be excluded when training CNV mutation rate models with `athena`. Curation steps described in `dsmap/annotations/curation/curate_athena_training_exclusion_mask.GRCh37.sh`. See below for the criteria for this file. |  
| Athena hg38 neutral region mask | `athena.hg38.training_exclusion.bed.gz` | _N/A_ | A list of regions under purifying selection to be excluded when training CNV mutation rate models with `athena`. Curation steps described in `dsmap/annotations/curation/curate_athena_training_exclusion_mask.hg38.sh`. See below for the criteria for this file. |  

---  

## Curation details  

### Athena hg38 neutral region mask  
This mask reflects the union of:  

1. Mammalian conserved regions per 91-way GERP (downloaded from Ensembl v103), and

2. Exons from genes with any marginal evidence of constraint against coding variants.  

Specifically, the exons in `2` were taken from any transcript in Gencode v37 with:  
*  observed/expected loss-of-function variants < 1 in gnomAD v2.1.1, or
*  observed/expected missense variants < 1 in gnomAD v2.1.1, or
*  probability of loss-of-function intolerance (pLI) > 0.1 in gnomAD v2.1.1

### Athena GRCh37 neutral region mask  
This mask was compiled identically to the hg38 mask (see above), except for using:    

1. 37-way GERP (downloaded from Ensembl v75), and

2. Gencode v19 (the last native GRCh37 Gencode annotation).  

