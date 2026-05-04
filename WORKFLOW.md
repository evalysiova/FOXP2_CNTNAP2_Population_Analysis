## **WORKFLOW - Population Genetics Analysis**

The full pipeline for Population Genetics Analysis using the 1000 Genomes Project Data.

### Requirements:
- **R**
- **R packages:** vcfR, adegenet, ade4, ggplot2, plotly, dplyr, tidyr, tidyverse, hierfstat, reshape2, pegas, GWLD, randomcoloR, stringdist, ape, 

-------------------------------------------------------------------------------

### Step 1 - SNPs Analysis in R

Population structure analysis and F-statistics of SNPs.

<br>

**Script availability:**
```r
source("scripts/snps_analysis.R")
```

<br>

> **Overview:**
>  - Loads data (all-snps.vcf.gz) in the environment through vcfR.
>  - Checks and produces unique variants IDs.
>  - Creates genind and genpop objects.
>  - Performs PCA, CA and DAPC.
>  - Visualizes contribution of contributing SNPs and their respective allele frequencies across 26 populations.
>  - Calculates Fst, Ho and He of the given genind object.
>  - Creates genind object for the contributing SNPs.
>  - Calculates Ho, He and per population/contributing SNPs +++
>  - Calculates Fis and Heterozygosity deviation per population (contributing variants).
>  - Performs population pairwise Fst and mean Fst-Ho Spearman correlation for the contributing SNPs.

-------------------------------------------------------------------------------

### Step 2 - INDELs Analysis in R

Population structure analysis and F-statistics of INDELs.

<br>

**Script availability:**
```r
source("scripts/indels_analysis.R")
```

<br>

> **Overview:**
>  - Loads data (all-indels.vcf.gz) in the environment through vcfR.
>  - Checks and produces unique variants IDs.
>  - Creates genind and genpop objects.
>  - Performs PCA, CA and DAPC.
>  - Visualizes contribution of contributing INDELs and their respective allele frequencies across 26 populations.
>  - Calculates Fst, Ho and He of the given genind object.
>  - Creates genind object for the contributing INDELs.
>  - Calculates Ho, He and per population/contributing INDEL +++
>  - Calculates Fis and Heterozygosity deviation per population (contributing variants).
>  - Performs population pairwise Fst and mean Fst-Ho Spearman correlation for the contributing INDELs.

-------------------------------------------------------------------------------

### Step 3 - Haplotypes Analysis in R

Population structure analysis and F-statistics of haplotypes consisting of the contributing variants. The dataset is handled as three distinct genomic loci,  one corresponding to FOXP2 coordinates and the two remaining to CNTNAP2 ones.

<br>

**Script availability:**
```r
source("scripts/haplotypes_analysis.R")
```

<br>

> **Overview:**
>  - Extraction of phased genotypes.
>  - Construction of haplotypes through phased data, handled as three loci.
>  - Identification of unique haplotypic sequence and calculation of frequency distribution per population. 
>  - Creation of a total genind object and three locus-specific ones.
>  - Calculation of Heterozygosity measures and performance of pairwise Fst (across the 26 populations) per locus.
>  - Sperman correlation of Ho and Fst per locus.
>  - Identifies population structure (PCA, CA, DAPC) and haplotypes contributing to population differentiation.
>  - Examines frequency distribution of haplotypes in locus 1, along with recoding of their respective sequence.
>  - Calculates hamming distances across the contributing haplotypes and performs PCoA.
>  - Extracts genotypes from the VCF files and converts them.
>  - Subsets per genomic locus.
>  - Calculates Linkage Disequilibrium of contributing variants (through Random Mutual Information) per population and visualizes them in LD plots.

-------------------------------------------------------------------------------

### Step 4 - Linkage Disequilibrium patterns around rs7782412 & rs7799652 

Exploration of LD patterns around the chromosomal region (+-50kb) that rs7782412 & rs7799652 SNPs reside. The used dataset includes SNPs and INDELs as distinct entries & genotypes are phased.
Filtering for common variants (universal MAF > 0.05) is performed as preprocess in the VCF file.

<br>

**Script availability:**
```r
source("scripts/LD_trait_locus.R")
```

<br>

> **Overview:**
> - Extraction of phased genotypes.
> - Calculation of RMI measure (GWLD package).
> - Visualization of patterns for the summary of variants, along with zoomed in heatmaps (20 variants).

-------------------------------------------------------------------------------

### ***NOTES***
 
 - SNPs and INDELs are analysed separately in the present pipeline (Steps 1-2).
 - All 26 populations are included in the present workflow.
 - Processed VCF files (all-snps.vcf.gz, etc.) are available in the present repository.
 - Preprocessing of the raw 1000 Genomes Project data was performed through BCFtools and Beagle 5.5. Reference scripts are available in the repository for replication.
