#/bin/bash

set -x


## Download 1000 Genomes data ##
# This file from the 1000 Genomes contains genetic data of 629 individuals from different ethnic backgrounds.
# Note, this file is quite large (>60 gigabyte).  
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz

# Convert vcf to Plink format.
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes --noweb
# Noteworthy, the file 'ALL.2of4intersection.20100804.genotypes.bim' contains SNPs without an rs-identifier, these SNPs are indicated with ".". This can also be observed in the file 'ALL.2of4intersection.20100804.genotypes.vcf.gz'. To check this file use this command: zmore ALL.2of4intersection.20100804.genotypes.vcf.gz .
# The missing rs-identifiers in the 1000 Genomes data are not a problem for this tutorial.
# However, for good practice, we will assign unique indentifiers to the SNPs with a missing rs-identifier (i.e., the SNPs with ".").
plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs --noweb

## QC on 1000 Genomes data.
# Remove variants based on missing genotype data.
plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS --noweb

# Remove individuals based on missing genotype data.
plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2 --noweb

# Remove variants based on missing genotype data.
plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3 --noweb

# Remove individuals based on missing genotype data.
plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS5 --noweb # Renaming output file to skip removing low frequency and rare MAF variants


