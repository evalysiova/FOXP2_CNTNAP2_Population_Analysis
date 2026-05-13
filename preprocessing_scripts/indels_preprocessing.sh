# Download the raw data from 1000 Genomes Project (GRCh37 version).
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# Keep regions of interest (FOXP2 & CNTNAP2 coordinates - GRCh37)

bcftools view -r 7:113721382-114338820,7:145808893-148123090 ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o genes_all.vcf.gz

# Keep only indels

bcftools view -v indels genes_all.vcf.gz -Oz -o indels.vcf.gz

# Please convert the files manually to GRCh38 (eg. Assembly Converter, Ensembl).
# Converted to GRCh38 => snps-lifted.vcf

# Extract variants (without IDs).

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' indels-lifted.vcf > variants.txt

# Annotate variants with IDs. The metadata were obtained from NCBI and an processed file is available in the repository.

bgzip matched_lines2.2.txt
tabix -s1 -b2 -e2 matched_lines2.2.txt.gz
bcftools annotate -a matched_lines2.2.txt.gz -c CHROM,POS,ID -o all-indels.vcf -Oz indels-lifted.vcf

# Check variants of given dataset.

bcftools view -H all-indels.vcf | wc -l     #3569 indels
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' all-indels.vcf > variants2.txt

# Check for multiallelic variants.

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' all-indels.vcf | awk -F',' '{ if (NF>1) print }' > multiallelic.txt
wc -l multiallelic.txt #210
