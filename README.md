## *FOXP2* & *CNTNAP2* Population Genetic Analysis

The present repository focuses on investigating single nucleotide (SNVs) and insertion/deletion variants (INDELs)
through the scope of Population Genetics, using 1000 Genomes Project Phase 3 data. All 26 populations and 5 superpopulations 
(AFR, EUR, EAS, SAS & AMR) available in the dataset are utilized.

-----------------------------------------------------
### Requirements
- R
- R packages:

### Optional Bioinformatics tools:
- BCFtools
- Beagle 5.5

------------------------------------------------------

### Data

The present analysis is conducted using the 1000 Genomes Project Data (Phase 3, GRCh37).
Original data are availble at:
https://www.internationalgenome.org/data-portal/data-collection/phase3

> Raw files are not available in this repository. VCF preprocessing can be regenerated through the available
> scripts (using Bash, BCFtools, Beagle 5.5 & Ensembl Assembly Converter). The analysis can be reproduced, additionally,
> through the publicly available files.

----------------------------------------------------------

### Repository Structure

```
├── README.md
├── WORKFLOW.md
├── PREPROCESSING_WORKFLOW.md
├── scripts/ 
│   ├── snps_analysis.R
│   ├── indels_analysis.R
│   ├── haplotypes_analysis.r
│   └── LD_trait_locus.R
├── preprocessing_scripts/
│   ├──
├── data
│   └── raw
├── output_snps
│   └── figures
├── output_indels
│   └── figures
└── output_LD
    └── figures

```
---------------------------------------------------------------
### Reproduce

#### - VCF preprocessing
The instructions for the according scripts are availble at PREPROCESSING_WORKFLOW.md.

#### - R Analysis
All approaches are listed in order and thoroughly explained at WORKFLOW.md of this repository.

All packages versions are recorder in `renv.lock`. The exact environment and version of packages can be restored with:

```r
install.packages("renv")
renv::restore()
```

> Processed VCF files can be directly downloaded through available commands in the R environment.
> The data can be found locally at generated directory ``data/raw/``.
> The total dataset constitutes of:
> - ``all-snps.vcf.gz``
> - ``all-indels.vcf.gz``
> - ``common_merged.vcf.gz``
> - ``contrib_ids.txt``
> - ``contrib_indels.vcf.gz``
> - ``contib_snps.vcf.gz``
> - ``contr_phased.vcf.gz``
> - ``pop.csv``

-----------------------------------------------------------------
### Output

All figures and files are saved through the working R environment.
Generated directories include:
- ``output_snps/``
- ``output_indels/``
- ``output_LD/``

-------------------------------------------------------------------
### Author

Evangelia-Sotiria Lysiova   (
evanlysi@mbg.duth.gr | 
evalysiova@gmail.com )
