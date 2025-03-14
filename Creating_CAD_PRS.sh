#################################################################################################
####################Script for Creating CAD Polygenic Risk Score in ALSPAC#######################
##SNPs here identified from summary data provided at https://doi.org/10.1038/s41588-022-01233-6##
#################################################################################################

#!/bin/bash

module unload -f compilers mpi gcc-libs
module load beta-modules
module load bgen/1.1.4
module load plink/2.0alpha-git
module load stata/15

#### Create bgen index files for each chromosome
##This bit takes quite a long time so should be submitted as a job (the rest can be run on log-in nodes if preferred)

for i in {01..22}; do bgenix -g /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/filtered_${i}.bgen -index; done

#### Now extract specific SNPs required for PRS 
##This step needs a .txt file containing a single column of the necessary SNP IDs
## ALSPAC bgen files don't have rsIDs and when using chr:pos bgenix seems to only be able to look for a range (e.g. XX:XXXX-XXXXX) rather than a single position
## Text file therefore altered before use here so that each SNP is in form of something like 01:123123-123123 for single SNP detection
## Also, this file seems to want a leading 0 before the 1-9 chromosomes so I had to add these too, but later plink wanted it removed again. Not sure why different.

for i in {01..22}; do bgenix -g /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/filtered_${i}.bgen -list -incl-range /home/rmgpstc/Scratch/ALSPAC/SNPs.txt; done
for i in {01..22}; do bgenix -g filtered_${i}.bgen -incl-range /home/rmgpstc/Scratch/ALSPAC/SNPs.txt > chr_${i}.bgen; done

#### Move new .bgen files containing only SNPs of interest to a different folder and concatenate into a single .bgen file containing all necessary SNPs for PRS
## This identifies 225 SNPs out of the 241 listed in summary stats 

mv chr_*.bgen /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/Converted_data/
cat-bgen -g chr_{01..22}.bgen -og /home/rmgpstc/Scratch/ALSPAC/Genotyping_data/Converted_data/PRS_SNPS.bgen

#### Calculate MAFs for SNPS in ALSPAC data to allow data harmonisation with GWAS data

plink2 --bgen PRS_SNPS.bgen --sample sample.sample --freq --out PRS_SNPS_maf
head PRS_SNPS_maf.afreq

#### Now perform allele harmonisation using Stata do file written for this purpose
##This will output a txt file called "Cleaned_GWAS_SNPs.txt" that can be used to generate the score in Plink

stata-xp
stata -b do "Harmonising_ALSPAC_Data.do"

#### Use PLINK to create PRS

## PLINK needs the score file to be saved in a space- or tab-delimited format or it returns an error.
## In the ALSPAC .bgen files the reference allele is listed as being first so have included ref-first option.
## The sample file supplied by ALSPAC also has to be modified to remove the mother and father columns before use as otherwise Plink gets confused (original swapped.sample -> sample.sample used here).
## The list-variants command tells you which variants were included in case some are missing. The ignore-dups removes duplicate SNPs with more than one alt allele.
## The cols=+scoresums option creates the actual PRS.
## 221 variants ultimately used out of 241 identified in the GWAS summary stats. Two were duplicates, two were palindromic (EAF 0.51) and so removed, and the others weren't present in the target data.

cd Converted_data/
plink2 --bgen PRS_SNPS.bgen ref-first --sample sample.sample --score Cleaned_GWAS_SNPs.txt list-variants ignore-dup-ids cols=+scoresums --out CAD_PRS