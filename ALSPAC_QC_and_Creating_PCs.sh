#################################################################################################
###########Script for QC of ALSPAC Genetic Data and Generation of Prinicpal Components###########
####################################Author: Scott Chiesa#########################################
#################################################################################################


#!/bin/bash


#### Load the various modules that may be needed for the script to run

module unload -f compilers mpi gcc-libs
module load beta-modules
module load plink/2.0alpha-git  ##This version of Plink on Myriad is old and can't handle some commands so I've copied the new version from GitHub into the folder of interest and am using that instead
module load stata/15


#===============Create one large file containing all autosomal chromosomes===============#

#### Convert all chromosome bgen files to pgen files

for i in {01..22}; do
  plink2 --bgen /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/filtered_${i}.bgen ref-first \
         --sample /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/Converted_data/sample.sample \
         --make-pgen \
         --out /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/PCA/chr_${i}
done

#### Change into directory to be used from now on

cd /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/PCA/

#### Remove multi-allelic variants as can't combine files with these present
##Identify multiallelic sites and print list of these to a text file for each chromosome

for i in {01..22}; do 
  awk '{print $3}' chr_${i}.pvar | sort | uniq -d > chr_${i}_dup.txt
done

## Now recreate pgen files but with these sites excluded

for i in {01..22}; do
 ./plink2 --pfile chr_${i} \
          --exclude chr_${i}_dup.txt \
          --make-pgen \
          --out chr_${i}_clean
done

#### Now merge these new files together into one big file

./plink2  --pfile ./chr_01_clean \
          --pmerge-list ./merge_list.txt \
          --make-pgen \
          --out ./merged_genome     ## txt file needed for this consists of single column containing chr_02_clean - chr_22_clean


#============Remove Unwanted People===============#

#### Remove mothers from sample as not needed

./plink2  --pfile ./merged_genome \
          --remove ./remove_mums.txt \
          --make-pgen \
          --out ./merged_genomes_nomums   ## remove_mums.txt consists of two columns containing FID and IID of mothers only (identified with ID ending F3M)

#### Remove siblings and non-white ethnicity

stata -b do "./Identifying_sibs_ethn.do"    ## this uses Stata to open the phenotype file and drop anyone with qlet == B (2nd child) and c800 != 1 (non-white) and then generates a txt file list of IDs for the remaining white non-related individuals

./plink2  --pfile ./merged_genomes_nomums \
          --keep ./remove_sib_ethn.txt \
          --make-pgen \
          --out ./merged_genomes_nomums_nosibs_white     ## txt file needed is that created in step above


#============Quality Control and Data Preparation===============#

### Most of this not necessary as ALSPAC data QCed before release but running anyway

#### Remove SNPs with low call rate (5%)

./plink2  --pfile ./merged_genomes_nomums_nosibs_white \
          --geno 0.05 \
          --make-pgen \
          --out ./merged_genomes_nomums_nosibs_white_lcr1

#### Remove participant samples with low call rate (5%)

./plink2  --pfile ./merged_genomes_nomums_nosibs_white_lcr1 \
          --mind 0.05 \
          --make-pgen \
          --out ./merged_genomes_nomums_nosibs_white_lcr2

#### Remove related individuals (< 3rd degree)

./plink2  --pfile ./merged_genomes_nomums_nosibs_white_lcr2 \
          --king-cutoff 0.125 \
          --make-pgen \
          --out ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated

#### Remove SNPs deviating from HWE (5x10e-07)

./plink2  --pfile ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated \
          --hwe 5e-7 \
          --make-pgen \
          --out ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe

#### Remove SNPs with MAF < 5%
./plink2  --pfile ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe \
          --maf 0.05 \
          --make-pgen \
          --out ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf


#============Prune Dataset and Check for Heterozygosity===============#

#### Create list of pruned SNPs not in LD
## Indep-pairwise needs a required window size in variant count, an optional variant count to shift the window at the end of each step, and a required r2 threshold. 
## At each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, and variants are greedily pruned from the window until no such pairs remain.

## Create lists of SNPs to be kept (prune.in) or discarded (prune.out)

./plink2  --pfile ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf \
          --indep-pairwise 250 50 0.1 \
          --out pruning_list

#### Remove participant samples with outlying heterozygosity (>3SDs from mean)
## Create het file that shows number of homozygous genotypes per person

./plink2  --pfile merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf \
          --extract pruning_list.prune.in \
          --het \
          --out het_check

 ## Now use Stata to identify any samples that demonstrate outlying heterozygosity
 
 stata -b do "./Identifying_outlying_het.do"    ## This uses Stata to calculate heterozygosity rates >3SDs from mean and identifies 26 people

 ## Now use PLINK to remove these people

./plink2  --pfile merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf \
          --remove het_outliers.txt \
          --make-pgen \
          --out merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet


#============Check Numbers at Every QC Stage===========#

wc -l merged_genomes.psam   ## n=17451
wc -l merged_genomes_nomums.psam    ##n=8798
wc -l merged_genomes_nomums_nosibs_white.psam   ## n=7795
wc -l merged_genomes_nomums_nosibs_white_lcr1.psam    ## n=7795
wc -l merged_genomes_nomums_nosibs_white_lcr2.psam    ## n=7795
wc -l merged_genomes_nomums_nosibs_white_lcr2_nonrelated.psam   ## n=7753
wc -l merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe.psam   ## n=7753
wc -l merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf.psam   ## n=7753
wc -l merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet.psam   ## n=7727


#============Generate Principal Components===========#

## Two files will be formed: .eigenval: eigenvalues explained by each of 10 PCs & .eigenvec: PCs for each individual

./plink2  --pfile ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet \
          --extract ./pruning_list.prune.in \
          --pca 10 \
          --out ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet_pruned

##### Check files to make sure all ok
## eigenval 
cat ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet_pruned.eigenval
wc -l ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet_pruned.eigenval  ##n=7727

## eigenvec
head ./merged_genomes_nomums_nosibs_white_lcr2_nonrelated_hwe_maf_nohet_pruned.eigenvec
