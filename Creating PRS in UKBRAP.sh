#################################################################################################
##################################PRS Construction in UKB RAP ###################################
####################################Author: Scott Chiesa 2025####################################
#################################################################################################

# Create new directory within session

mkdir ./tools

# Download required programmes 

wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20250129.zip #PLINK2
wget https://www.chg.ox.ac.uk/~gav/resources/bgen_v1.1.4-Ubuntu16.04-x86_64.tgz #Bgenix
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip #PRSice

# Unzip all to tools folder and give access rights where necessary

unzip ./plink2_linux_* -d ./tools
chmod a+x ./tools/plink2
tar -xzf ./bgen_v1.1.4-Ubuntu16.04-x86_64.tgz --strip-components=1 -C ./tools
unzip ./PRSice_linux -d ./tools


#################################################################################################
#####################################Preparing Target Data ######################################
#################################################################################################


#===============Initial QC 1 on individual BGEN and then conversion to PGEN==============#

  # Remove SNPs with low call rate (5%)
  # Remove SNPs deviating from HWE (5x10e-07)
  # Remove SNPs with MAF <1%

for i in {1..22}; do 
./tools/plink2 --bgen "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${i}_b0_v3.bgen" ref-first \
               --sample "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${i}_b0_v3.sample" \
               --geno 0.05 \
               --hwe 5e-7 \
               --maf 0.01 \
               --make-pgen \
               --out ./chr_${i}_qc1


#=======Remove multi-allelic variants as can't combine files with these present==========#

##Identify multiallelic sites and print list of these to a text file for each chromosome

for i in {1..22}; do 
  awk '{print $3}' chr_${i}_qc1.pvar | sort | uniq -d > chr_${i}_qc1_dup.txt
done

## Now recreate pgen files again but with these sites excluded

for i in {1..22}; do
./tools/plink2 	--pfile ./chr_${i}_qc1 \
        	  	  --exclude ./chr_${i}_qc1_dup.txt \
          		  --make-pgen \
          		  --out ./chr_${i}_qc1_clean
done

#### Now merge these new files together into one big file

./tools/plink2 	--pfile ./chr_1_qc1_clean \
          		  --pmerge-list /mnt/project/Scripts/merge_list.txt \
          		  --make-pgen \
          		  --out ./merged_qc1_clean     ## the txt file needed above simply consists of single column containing chr_2_qc_clean - chr_22_qc_clean

## Check numbers 

wc -l ./merged_qc1_clean.psam   ## n=

#===============Data QC 2 on Sample==============#

# Remove participant samples with low call rate (5%)
# Identify and remove palindromic SNPs

  ## Remove samples with low call rate
	
./tools/plink2 --pfile ./merged_qc1_clean \
               --freq \
  				     --mind 0.05 \
               --make-pgen \
  				     --out ./merged_qc2_clean

   ## Check numbers 

wc -l ./merged_qc2_clean.psam   ## n=             


#===============Data QC 3 on Sample==============#


  ## Identify and remove palindromic SNPs with allele frequencies close to 0.5 (0.45 - 0.55)

awk '/^[^#]/ { if( $5>0.45 && $5<0.55 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' \
  merged_qc2_clean.afreq > excl_ambiguous.txt

./tools/plink2 --pfile ./merged_qc2_clean \
               --remove ./excl_ambiguous.txt 
               --mind 0.05 \
               --make-pgen \
               --out ./merged_qc3_clean

 ## Check numbers 

wc -l ./merged_qc3_clean.psam   ## n=

#===============Data QC 4 Before Creating PRS==============#

  # Removing participants with non-white ethnicity, relatives, or sex discordance
  # This can easily be modified if you want to include everyone

    ## First create a data dictionary with all available phenotypes

dx extract_dataset project-J1qvq8QJ69pGx77pKyKFb915:record-J1qyKFQJjX0yK25BQ777YXgF -ddd --delimiter ","

    ## Then extract the necessary fields for QC
    ## These are ID, reported sex, genetic sex, white ancestry, chromosomal aneuploidy, and whether they were used by UKB to calculate PCs (which means they're not related to anyone else)

dx extract_dataset project-J1qvq8QJ69pGx77pKyKFb915:record-J1qyKFQJjX0yK25BQ777YXgF \
    --fields participant.eid,participant.p31,participant.p22001,participant.p22006,participant.p22019,participant.p22020 \
    -o ./qc4_vars.csv

# Now create a file to merge with genotypes while filtering out undesirables
    ## Final output is csv with FID, IID, sex, genetic sex, white ancestry, no aneuploidy, no relatives, and all other phenotypes

# Run R

R

# Set paths to R packages

.libPaths(c("/home/dnanexus/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

# Load library needed

library(dplyr)

# Read in the csv generated above

df <- read.csv("./qc4_vars.csv", stringsAsFactors = FALSE)

# Create flags for failures
    ## This checks if reported sex matches genetic, if ancestry is recorded as White-British (1), if aneuploidy is reported, and if they were included in generation of PCs

df <- df %>%
  mutate(
    sex_mismatch = participant.p31 != participant.p22001,
    non_white_british = participant.p22006 != 1 | is.na(participant.p22006),
    aneuploidy = !is.na(participant.p22019),
    related = participant.p22020 != 1
  )

# Count how many fail each criterion and print list
    ## Lets us see why people were removed

qc_summary <- df %>%
  summarise(
    total = n(),
    fail_sex_mismatch = sum(sex_mismatch, na.rm = TRUE),
    fail_non_white_british = sum(non_white_british, na.rm = TRUE),
    fail_aneuploidy = sum(aneuploidy, na.rm = TRUE),
    fail_related = sum(related, na.rm = TRUE),
    pass_all = sum(!sex_mismatch & !non_white_british & !aneuploidy & !related, na.rm = TRUE)
  )

print(qc_summary)

# Now apply QC filters for necessary output file

df_qced <- df %>%
  filter(
    participant.p31 == participant.p22001, # reported sex == genetic sex
    participant.p22006 == 1,               # White British ancestry
    is.na(participant.p22019),             # No sex chromosome aneuploidy
    participant.p22020 == 1                # Used in PCA (unrelated)
  )

# Rename columns

df_qced <- df_qced %>%
  rename(
    IID = participant.eid,
    sex = participant.p31,
    genetic_sex = participant.p22001,
    white_ancestry = participant.p22006,
    no_aneuploidy = participant.p22019,
    no_relatives = participant.p22020
  )

# Add FID column (same as IID)
df_qced$FID <- df_qced$IID

# Select columns for phenotype table
cols <- c("FID", "IID", "sex", "genetic_sex", "white_ancestry", "no_aneuploidy", "no_relatives")
df_phenotype <- df_qced[, cols]

# Save output
write.csv(df_phenotype, "./phenotypes_for_qc4.csv", row.names = FALSE)

# Quit R

q()


#============Create new pgen file with only these people identified above===========#

# Take csv file and keep only IID and FID with no headers to use in plink command

awk -F, 'NR>1 {print $1, $2}' phenotypes_for_qc4.csv > qced_keep.txt

# Use plink to recreate final merged QCed pgen dataset that contain only these people

.tools/plink2 --pfile merged_qc3_clean \
  --keep qced_keep.txt \
  --make-pgen \
  --out final_qced_dataset


#################################################################################################
######################################Preparing Base Data #######################################
#################################################################################################


#============Prepare GWAS File to be used in PRSice=================#

## Download Summary Stats (in this case Evangelou et al SBP stats but can adapt to do multiple)

wget -O SBP.txt.gz https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz

## Run R

R

## Set paths to saved R packages in case needed

.libPaths(c("/home/dnanexus/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

##Load library needed

library(data.table)

## Read in .gz file, rename columns, drop those not needed, save as tab-delimited txt file

gwas_data <- fread("./SBP.txt.gz", header=T) ## now check columns using 'head(gwas_data, 100)'
names(gwas_data) <- c("CHRBP","A1","A2","EAF","BETA","SE","P","cases","effective_cases")
gwas_data[, c("cases","effective_cases") := NULL]
fwrite(gwas_data, "./sbp_gwas_data.txt", sep = " ")

## Quit R

q()

## Modify columns in txt file so ready for PRSice
    # This will modify the first column (e.g. 10:123456:SNP) into three columns called SNP / CHR / and BP (e.g. 10:123456 / 10 / and 123456)
    # The rest stay as is


awk '
BEGIN {OFS="\t"}
NR==1 {
  print "SNP","CHR","BP",$2,$3,$4,$5,$6,$7
  next
}
{
  split($1,a,":")        
  snp=a[1]":"a[2]        
  chr=a[1]              
  pos=a[2]               
  print snp,chr,pos,$2,$3,$4,$5,$6,$7
}
' sbp_gwas_data.txt > sbp_gwas_final.txt


#################################################################################################
##########################################Creating PRS###########################################
#################################################################################################


#=======================Use PRSice to create scores===========================#

Rscript ./tools/PRSice.R \
  --prsice ./tools/PRSice_linux \
  --base ./sbp_gwas_final.txt \
  --target ./final_qced_dataset \
  --stat BETA \
  --snp SNP \
  --chr CHR \
  --bp BP \
  --A1 A1 \
  --A2 A2 \
  --pvalue P \
  --score std \
  --bar-levels 5e-8,1e-6,1e-4,1e-2,0.05,0.1,0.5,1 \
  --fastscore \
  --out ./SBP_PRS \
  --clump-kb 250 \
  --clump-r2 0.1 \
  --clump-p 1.0 \
  --extract ./opt/notebooks/PRS/trait_prs.valid \
  --no-regress \
  --binary-target T          