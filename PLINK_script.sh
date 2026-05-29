#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --mem=8G

module load plink/1.90b6.21
module load r/3.6.0


#converted output to the file types needed. 
#─plink --file /work/users/t/y/tycook/Jackson_ancestry/GWA/GitHub/PLINK_091024_1033/Jackson_01 --make-bed --out /work/users/t/y/tycook/Jackson_ancestry/GWA/attmpt_2/Jackson_01_binary
#above genereated an error about invalid chromosome so i did the below. 
plink --file /work/users/t/y/tycook/Jackson_ancestry/GWA/GitHub/PLINK_091024_1033/Jackson_01 --make-bed --out /work/users/t/y/tycook/new_ancestry_rerun/Jackson_01_binary --allow-extra-chr

# Investigate missingness per individual and per SNP and make histograms.
plink --bfile /work/users/t/y/tycook/Jackson_ancestry/GWA/Jackson_01_binary --missing --allow-extra-chr
# output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.


# Generate plots to visualize the missingness results.
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/hist_miss.R

# Delete SNPs and individuals with high levels of missingness, explanation of this and all following steps can be found in box 1 and table 1 of the article mentioned in the comments of this script.
# The following two QC commands will not remove any SNPs or individuals. However, it is good practice to start the QC with these non-stringent thresholds.  
# Delete SNPs with missingness >0.2.
plink --bfile Jackson_01_binary --geno 0.2 --make-bed --out Jackson_01_binary_rmv_snps_20 --allow-extra-chr

# Delete individuals with missingness >0.2.
plink --bfile Jackson_01_binary_rmv_snps_20 --mind 0.2 --make-bed --out Jackson_01_binary_missingness_20 --allow-extra-chr

# Delete SNPs with missingness >0.02.
plink --bfile Jackson_01_binary_missingness_20 --geno 0.02 --make-bed --out Jackson_01_binary_snps_missingness_2 --allow-extra-chr

# Delete individuals with missingness >0.02.
plink --bfile Jackson_01_binary_snps_missingness_2 --mind 0.02 --make-bed --out Jackson_01_binary_missingness_2 --allow-extra-chr





###################################################################
### Step2 ####

# Check for sex discrepancy.
# Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

plink --bfile Jackson_01_binary_missingness_2 --check-sex 

# Generate plots to visualize the sex-check results.
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/gender_check.R

# These checks indicate that there is one woman with a sex discrepancy, F value of 0.99. (When using other datasets often a few discrepancies will be found). 

# The following two scripts can be used to deal with individuals with a sex discrepancy.
# Note, please use one of the two options below to generate the bfile hapmap_r23a_6, this file we will use in the next step of this tutorial.

# 1) Delete individuals with sex discrepancy.
#grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
# This command generates a list of individuals with the status �PROBLEM�.
plink --bfile Jackson_01_binary_missingness_2 --make-bed --out Jackson_01_binary_missingness_2_sexdisc
# This command removes the list of individuals with the status �PROBLEM�.

# 2) impute-sex.
#plink --bfile HapMap_3_r3_5 --impute-sex --make-bed --out HapMap_3_r3_6
# This imputes the sex based on the genotype information into your data set.

###################################################
### Step 3 ### 

# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' Jackson_01_binary_missingness_2_sexdisc.bim > snp_1_22.txt
plink --bfile Jackson_01_binary_missingness_2_sexdisc --extract snp_1_22.txt --make-bed --out Jackson_01_binary_missingness_2_sexdisc_autosomal


# Generate a plot of the MAF distribution.
plink --bfile Jackson_01_binary_missingness_2_sexdisc_autosomal --freq --out MAF_check
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/MAF_check.R

# Remove SNPs with a low MAF frequency.
plink --bfile Jackson_01_binary_missingness_2_sexdisc_autosomal --maf 0.05 --make-bed --out Jackson_01_binary_missingness_2_4_low_MAF_frequency
# 1073226 SNPs are left
# A conventional MAF threshold for a regular GWAS is between 0.01 or 0.05, depending on sample size.


####################################################
### Step 4 ###

# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
# Check the distribution of HWE p-values of all SNPs.

plink --bfile Jackson_01_binary_missingness_2_4_low_MAF_frequency --hardy
# Selecting SNPs with HWE p-value below 0.00001, required for one of the two plot generated by the next Rscript, allows to zoom in on strongly deviating SNPs. 
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/hwe.R

# By default the --hwe option in plink only filters for controls.
# Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data.
plink --bfile Jackson_01_binary_missingness_2_4_low_MAF_frequency --hwe 1e-6 --make-bed --out Jackson_hwe_filter_step1

# The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE. 
# This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed
plink --bfile Jackson_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out Jackson_hwe_filter_step2

# Theoretical background for this step is given in our accompanying article: https://www.ncbi.nlm.nih.gov/pubmed/29484742 .

############################################################
### step 5 ###

# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions (inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise�.
# The parameters �50 5 0.2� stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.


plink --bfile Jackson_hwe_filter_step2 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
# Note, don't delete the file indepSNP.prune.in, we will use this file in later steps of the tutorial.
plink --bfile Jackson_hwe_filter_step2 --extract indepSNP.prune.in --het --out R_check
#### This file contains your pruned data set.
did not run
##Plot of the heterozygosity rate distribution
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/check_heterozygosity_rate.R

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# For data manipulation we recommend using UNIX. However, when performing statistical calculations R might be more convenient, hence the use of the Rscript for this step:
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/heterozygosity_outliers_list.R


# Output of the command above: fail-het-qc.txt .
# When using our example data/the HapMap data this list contains 2 individuals (i.e., two individuals have a heterozygosity rate deviating more than 3 SD's from the mean).
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.
plink --bfile Jackson_hwe_filter_step2 --remove het_fail_ind.txt --make-bed --out Jackson_hetero_outliers_rmd


############################################################
### step 6 ###

# It is essential to check datasets you analyse for cryptic relatedness.
# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial.

# Check for relationships between individuals with a pihat > 0.2.
plink --bfile Jackson_hetero_outliers_rmd --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

# The HapMap dataset is known to contain parent-offspring relations. 
# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome

# Generate a plot to assess the type of relationship.
Rscript --no-save /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/Relatedness.R

# The generated plots show a considerable amount of related individuals (explentation plot; PO = parent-offspring, UN = unrelated individuals) in the Hapmap data, this is expected since the dataset was constructed as such.
# Normally, family based data should be analyzed using specific family based methods. In this tutorial, for demonstrative purposes, we treat the relatedness as cryptic relatedness in a random population sample.
# In this tutorial, we aim to remove all 'relatedness' from our dataset.
# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).

plink --bfile Jackson_hetero_outliers_rmd --filter-founders --make-bed --out Jackson_hetero_founders

# Now we will look again for individuals with a pihat >0.2.
plink --bfile Jackson_hetero_founders --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
# The file 'pihat_min0.2_in_founders.genome' shows that, after exclusion of all non-founders, only 1 individual pair with a pihat greater than 0.2 remains in the HapMap data.
# This is likely to be a full sib or DZ twin pair based on the Z values. Noteworthy, they were not given the same family identity (FID) in the HapMap data.

# For each pair of 'related' individuals with a pihat > 0.2, we recommend to remove the individual with the lowest call rate. 
plink --bfile Jackson_hetero_founders --missing
# Use an UNIX text editor (e.g., vi(m) ) to check which individual has the highest call rate in the 'related pair'. 

# Generate a list of FID and IID of the individual(s) with a Pihat above 0.2, to check who had the lower call rate of the pair.
# In our dataset the individual 13291  NA07045 had the lower call rate.
vi 0.2_low_call_rate_pihat.txt
i 
13291  NA07045
# Press esc on keyboard!
:x
# Press enter on keyboard
# In case of multiple 'related' pairs, the list generated above can be extended using the same method as for our lone 'related' pair.

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
plink --bfile Jackson_hetero_founders --remove 0.2_low_call_rate_pihat.txt --make-bed --out Jackson_hetero_founders_low_callrate

################################################################################################################################

# CONGRATULATIONS!! You've just succesfully completed the first tutorial! You are now able to conduct a proper genetic QC. 


######      POPULATION STRATIFICATION!!!!! #####################

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


## Download 1000 Genomes data ##
# This file from the 1000 Genomes contains genetic data of 629 individuals from different ethnic backgrounds.
# Note, this file is quite large (>60 gigabyte).  
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz

# Convert vcf to Plink format.
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE

#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE

#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE#######STOP HERE
#######STOP HERE#######STOP HEREV#######STOP HERE#######STOP HERE#######STOP HERE
# Noteworthy, the file 'ALL.2of4intersection.20100804.genotypes.bim' contains SNPs without an rs-identifier, these SNPs are indicated with ".". This can also be observed in the file 'ALL.2of4intersection.20100804.genotypes.vcf.gz'. To check this file use this command: zmore ALL.2of4intersection.20100804.genotypes.vcf.gz .
# The missing rs-identifiers in the 1000 Genomes data are not a problem for this tutorial.
# However, for good practice, we will assign unique indentifiers to the SNPs with a missing rs-identifier (i.e., the SNPs with ".").
plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs

## QC on 1000 Genomes data.
# Remove variants based on missing genotype data.
plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS

# Remove individuals based on missing genotype data.
plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2

# Remove variants based on missing genotype data.
plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3

# Remove individuals based on missing genotype data.
plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4

# Remove variants based on MAF.
plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5

# Extract the variants present in HapMap dataset from the 1000 genomes dataset.
awk '{print$2}' Jackson_hetero_founders_low_callrate.bim > HapMap_SNPs.txt
plink --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6

# Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
plink --bfile Jackson_hetero_founders_low_callrate --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out Jackson_hetero_founders_low_callrate_MDS
# The datasets now contain the exact same variants.

## The datasets must have the same build. Change the build 1000 Genomes data build.
awk '{print$2,$4}' Jackson_hetero_founders_low_callrate_MDS.map > buildhapmap.txt
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
plink --bfile Jackson_hetero_founders_low_callrate_MDS --reference-allele 1kg_ref-list.txt --make-bed --out Jackson_hetero_founders_low_callrate_adj
# The 1kG_MDS7 and the HapMap-adj have the same reference genome for all SNPs.
# This command will generate some warnings for impossible A1 allele assignment.

# 2) Resolve strand issues.
# Check for potential strand issues.
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' Jackson_hetero_founders_low_callrate_adj.bim > Jackson_hetero_founders_low_callrate_adj_tmp
sort 1kGMDS7_tmp Jackson_hetero_founders_low_callrate_adj_tmp |uniq -u > all_differences.txt
# 1624 differences between the files, some of these might be due to strand issues.

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
# Generates a file of 812 SNPs. These are the non-corresponding SNPs between the two files. 
# Flip the 812 non-corresponding SNPs. 
plink --bfile Jackson_hetero_founders_low_callrate_adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap

# Check for SNPs which are still problematic after they have been flipped.
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
# This file demonstrates that there are 84 differences between the files.

# 3) Remove problematic SNPs from HapMap and 1000 Genomes.
#do noot have to do at present
#awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
# The command above generates a list of the 42 SNPs which caused the 84 differences between the HapMap and the 1000 Genomes data sets after flipping and setting of the reference genome.

# Remove the 42 problematic SNPs from both datasets.
#plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS2
#plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8

# Merge HapMap with 1000 Genomes Data.
plink --bfile corrected_hapmap --bmerge 1kG_MDS7.bed 1kG_MDS7.bim 1kG_MDS7.fam --allow-no-sex --make-bed --out MDS_merge2

# Note, we are fully aware of the sample overlap between the HapMap and 1000 Genomes datasets. However, for the purpose of this tutorial this is not important.

## Perform MDS on HapMap-CEU data anchored by 1000 Genomes data.
# Using a set of pruned SNPs
plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2
plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2
#plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --pca 10 --out pca_merge2

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
awk '{print$1,$2,"OWN"}' corrected_hapmap.fam>racefile_own.txt

# Concatenate racefiles.
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

# Generate population stratification plot.
Rscript /work/users/t/y/tycook/new_ancestry_rerun/plink_scripts/MDS_merged.R
# The output file MDS.pdf demonstrates that our �own� data falls within the European group of the 1000 genomes data. Therefore, we do not have to remove subjects.
# For educational purposes however, we give scripts below to filter out population stratification outliers. Please execute the script below in order to generate the appropriate files for the next tutorial.

## Exclude ethnic outliers.
# Select individuals in HapMap data below cut-off thresholds. The cut-off levels are not fixed thresholds but have to be determined based on the visualization of the first two dimensions. To exclude ethnic outliers, the thresholds need to be set around the cluster of population of interest.
#awk '{ if ($4 <-0.04 && $5 >0.03) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2

# Extract these individuals in HapMap data.
#plink --bfile HapMap_3_r3_12 --keep EUR_MDS_merge2 --make-bed --out HapMap_3_r3_13
# Note, since our HapMap data did include any ethnic outliers, no individuls were removed at this step. However, if our data would have included individuals outside of the thresholds we set, then these individuals would have been removed.

## Create covariates based on MDS.
# Perform an MDS ONLY on HapMap data without ethnic outliers. The values of the 10 MDS dimensions are subsequently used as covariates in the association analysis in the third tutorial.
#plink --bfile HapMap_3_r3_13 --extract indepSNP.prune.in --genome --out HapMap_3_r3_13
#plink --bfile HapMap_3_r3_13 --read-genome HapMap_3_r3_13.genome --cluster --mds-plot 10 --out HapMap_3_r3_13_mds

# Change the format of the .mds file into a plink covariate file.
#awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' HapMap_3_r3_13_mds.mds > covar_mds.txt

# The values in covar_mds.txt will be used as covariates, to adjust for remaining population stratification, in the third tutorial where we will perform a genome-wide association analysis.

##########################################################################################################################################################################

## CONGRATULATIONS you have succesfully controlled your data for population stratification!

# For the next tutorial you need the following files:
# - HapMap_3_r3_13 (the bfile, i.e., HapMap_3_r3_13.bed,HapMap_3_r3_13.bim,and HapMap_3_r3_13.fam
# - covar_mds.txt