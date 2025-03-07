#!/bin/bash
set -xu
################ Explanation of the main script ##########################

# This is the main script for second tutorial from our comprehensive tutorial on GWAS and PRS.
# To run this script the following (b)files from the first tutorial are required: HapMap_3_r3_12 (this bfile contain: HapMap_3_r3_12.fam,HapMap_3_r3_12.bim, and HapMap_3_r3_12.bed; you need all three), and indepSNP.prune.in.
# In this tutorial we are going to check for population stratification.
# We will do this as follows, the bfile (HapMap_3_r3_12) generated at the end of the previous tutorial (1_QC_GWAS) is going to checked for population stratification using data from the 1000 Genomes Project. Individuals with a non-European ethnic background will be removed.
# Furthermore, this tutorial will generate a covariate file which helps to adust for remaining population stratification within the European subjects.
# In order to complete this tutorial it is necessary to have generated the bfile 'HapMap_3_r3_12' and the file 'indepSNP.prune.in' from the previous tutorial.


##############################################################
############### START ANALISIS ###############################
##############################################################

# Copy the (b)files from the previous tutorial to the current directory (see explanation of the main script).
cp ../1_QC_GWAS/HapMap_3_r3_12.*  .
cp ../1_QC_GWAS/indepSNP.prune.in . 


## Download 1000 Genomes data ##
# This file from the 1000 Genomes contains genetic data of 629 individuals from different ethnic backgrounds.
# Note, this file is quite large (>60 gigabyte).  
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz

# Convert vcf to Plink format.
#plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
# Noteworthy, the file 'ALL.2of4intersection.20100804.genotypes.bim' contains SNPs without an rs-identifier, these SNPs are indicated with ".". This can also be observed in the file 'ALL.2of4intersection.20100804.genotypes.vcf.gz'. To check this file use this command: zmore ALL.2of4intersection.20100804.genotypes.vcf.gz .
# The missing rs-identifiers in the 1000 Genomes data are not a problem for this tutorial.
# However, for good practice, we will assign unique indentifiers to the SNPs with a missing rs-identifier (i.e., the SNPs with ".").
#plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs

## QC on 1000 Genomes data.
# Remove variants based on missing genotype data.
#plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS

# Remove individuals based on missing genotype data.
#plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2

# Remove variants based on missing genotype data.
#plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3

# Remove individuals based on missing genotype data.
#plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS5 # no maf filtering

# Remove variants based on MAF.
#plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5

# Extract the variants present in HapMap dataset from the 1000 genomes dataset.
awk '{print$2}' HapMap_3_r3_12.bim > HapMap_SNPs.txt
plink --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6

# Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
plink --bfile HapMap_3_r3_12 --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out HapMap_MDS
# The datasets now contain the exact same variants.

## The datasets must have the same build. Change the build 1000 Genomes data build.
awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
# buildhapmap.txt contains one SNP-id and physical position per line.

plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
# 1kG_MDS7 and HapMap_MDS now have the same build.

## Merge the HapMap and 1000 Genomes data sets

# Prior to merging 1000 Genomes data with the HapMap data we want to make sure that the files are mergeable, for this we conduct 3 steps:
# 1) Make sure the reference genome is similar in the HapMap and the 1000 Genomes Project datasets.
# 2) Resolve strand issues.
# 3) Remove the SNPs which after the previous two steps still differ between datasets.

# The following steps are maybe quite technical in terms of commands, but we just compare the two data sets and make sure they correspond.

# 1) set reference genome 
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
plink --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj
# The 1kG_MDS7 and the HapMap-adj have the same reference genome for all SNPs.
# This command will generate some warnings for impossible A1 allele assignment.

# 2) Resolve strand issues.
# Check for potential strand issues.
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp
sort 1kGMDS7_tmp HapMap-adj_tmp |uniq -u > all_differences.txt
# 1624 differences between the files, some of these might be due to strand issues.

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
# Generates a file of 812 SNPs. These are the non-corresponding SNPs between the two files. 
# Flip the 812 non-corresponding SNPs. 
plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap

# Check for SNPs which are still problematic after they have been flipped.
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
# This file demonstrates that there are 84 differences between the files.

# 3) Remove problematic SNPs from HapMap and 1000 Genomes.
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
# The command above generates a list of the 42 SNPs which caused the 84 differences between the HapMap and the 1000 Genomes data sets after flipping and setting of the reference genome.

# Remove the 42 problematic SNPs from both datasets.
plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS2
plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8

# Merge HapMap with 1000 Genomes Data.
plink --bfile HapMap_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2

# Note, we are fully aware of the sample overlap between the HapMap and 1000 Genomes datasets. However, for the purpose of this tutorial this is not important.

## Perform MDS on HapMap-CEU data anchored by 1000 Genomes data.
# Using a set of pruned SNPs
plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2
plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2


### MDS-plot

# Download the file with population information of the 1000 genomes dataset.
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel
# The file 20100804.ALL.panel contains population codes of the individuals of 1000 genomes.

# Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN, and EUR).
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt

# Create a racefile of your own data.
awk '{print$1,$2,"OWN"}' HapMap_MDS.fam>racefile_own.txt

# Concatenate racefiles.
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

# Generate population stratification plot.
Rscript MDS_merged.R 
# The output file MDS.pdf demonstrates that our �own� data falls within the European group of the 1000 genomes data. Therefore, we do not have to remove subjects.
# For educational purposes however, we give scripts below to filter out population stratification outliers. Please execute the script below in order to generate the appropriate files for the next tutorial.

## Exclude ethnic outliers.
# Select individuals in HapMap data below cut-off thresholds. The cut-off levels are not fixed thresholds but have to be determined based on the visualization of the first two dimensions. To exclude ethnic outliers, the thresholds need to be set around the cluster of population of interest.
awk '{ if ($4 <-0.04 && $5 >0.03) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2

# Extract these individuals in HapMap data.
plink --bfile HapMap_3_r3_12 --keep EUR_MDS_merge2 --make-bed --out HapMap_3_r3_13
# Note, since our HapMap data did include any ethnic outliers, no individuls were removed at this step. However, if our data would have included individuals outside of the thresholds we set, then these individuals would have been removed.

## Create covariates based on MDS.
# Perform an MDS ONLY on HapMap data without ethnic outliers. The values of the 10 MDS dimensions are subsequently used as covariates in the association analysis in the third tutorial.
plink --bfile HapMap_3_r3_13 --extract indepSNP.prune.in --genome --out HapMap_3_r3_13
plink --bfile HapMap_3_r3_13 --read-genome HapMap_3_r3_13.genome --cluster --mds-plot 10 --out HapMap_3_r3_13_mds

# Change the format of the .mds file into a plink covariate file.
awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' HapMap_3_r3_13_mds.mds > covar_mds.txt

# The values in covar_mds.txt will be used as covariates, to adjust for remaining population stratification, in the third tutorial where we will perform a genome-wide association analysis.

##########################################################################################################################################################################

## CONGRATULATIONS you have succesfully controlled your data for population stratification!

# For the next tutorial you need the following files:
# - HapMap_3_r3_13 (the bfile, i.e., HapMap_3_r3_13.bed,HapMap_3_r3_13.bim,and HapMap_3_r3_13.fam
# - covar_mds.txt





