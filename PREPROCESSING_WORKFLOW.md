## PREPROCESSING WORKFLOW

The full pipeline, for the raw 100 Genomes Project VCF (chromosome 7) to FOXP2 and CNTNAP2 specific datasets.

<md>
  
This workflow is optional for reproducibility. All processed VCF files are available in Releases of the 
present repository. It is recommended to run the scripts line by line, as some steps require visual inspections.
Additionally, liftover of data (GRCh37 to CRCh28) should be performed manually with the tool of choice. 
Herein Assembly Converter of Ensembl was used, but CrossMap or other online tools can be utilized.

### Requirements:
- Bash
- BCFtools
- Beagle 5.5

------------------------------------------------

### Step 0 - Download Raw Data
Download phased VCF files from 1000 Genomes Project (Phase3 - GRCh37):
https://www.internationalgenome.org/data/

<md>

> Data can be downloaded through the terminal with:
>  <md>
>
>`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
>
>`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi`
>
> <md>
>
> The commands are also availble in the `snps_preprocessing.txt` script.

-----------------------------------------------

### Step 1 - SNPs Preprocessing

Raw 1000 Genomes Project to SNV gene-specific dataset.

#### Script availability:
`snps_preprocessing.txt`

> Overview:
> - Subsets for the chromosomal coordinates of interest (FOXP2 and CNTNAP2 in GRCh37) within three partial files.
>  - Excludes INDELs within the dataset.
>  - Utilizes the converted VCF files (from Assembly Converter of Ensembl in 
    GRCh38) to merge the three parts.
>  - Annotates the variants with rs identifiers.
>  - Extracts a list of all present SNPs and their respective coordinates and alleles.
>  - Checks for multiallelic markers.

--------------------------------------

### Step 2 - INDELs preprocessing

Raw 1000 Genomes Project to INDELs gene-specific data.

#### Script availability:
`indels_preprocessing.txt`

> Overview:
>  - Subsets for the chromosomal coordinates of interest (FOXP2 and CNTNAP2).
>  - Exludes SNPs within the dataset.
>  - Utilizes the converted VCF file (from Assembly Converter of Ensembl in GRCh38) to annotate with identifiers.
>  - Extracts a list with the summary of the INDELs and their respective coordinates and alleles.
>  - Checks for multiallelic variants.

--------------------------------------

Step 3 - Contributing variants processing

Gene-specific to contributing variants-specific datasets.
> This step should be performed following the DAPC step of SNP & INDEL's Analysis in R (see WORKFLOW.md, step 1-2).

#### Script availability:
`contrib_vcfs_scripts.txt`

> Overview:
> - Subesets the all-snps.vcf & all-indels.vcf using the contrib_snps_sorted.txt & contrib_indels.txt listsas a guide.

--------------------------------------------

Step 4 - Haplotypes Preprocessing

Unphased variant-specific to phased data of contributing variants.

#### Script availability:
`haplotypes_vcf_preprocessing.txt`

> Overview:
> 
>  - Subsets snps-all.vcf and indels.vcf.gz for variants in chromosomal coordinates of interest (loci 1-3)
>  - Handles SNPs and INDELs with identical coordinates as distinct records.
>  - Phases data with BEAGLE 5.5 and Java.
>  - Adds INFO column in the phased files.
>  - Subsets for only the contributing variants.

----------------------------------------------
Step 5 - LD preprocessing

VCF files with all variants upstream and downstream of rs7782412 & rs7799652.

#### Script availability:
`trait_locus_vcf_preprocessing.txt`

> Overview:
> - Subsets for 50 kb upstream and downstream of the variants of interest.
> - Handles SNPs and INDELs with identical coordinates as distinct records.
> - Phases data with BEAGLE 5.5 and Java.
> - Adds INFO column in the phased files.
> - Filters for common variants (MAF > 0.05).

