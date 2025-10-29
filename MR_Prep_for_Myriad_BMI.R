# Performing a two-sample MR between education and BMI in UK Biobank using R package TwoSampleMR #
# Work in progress - needs updating#

## The below packages need to be installed in the first instance so that the subsequent libraries can be loaded.

#install.packages("remotes")
#install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
#install.packages("genetics.binaRies", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
#install.packages(c("tidyverse", "data.table", "meta", "devtools", "pacman", "ieugwasr", "LDlinkR", "MRPRESSO", "extrafont"))
#install.packages("devtools")
#remotes::install_github("phenoscanner/phenoscanner")
#devtools::install_github("MRCIEU/MRInstruments")

.libPaths(c("~/ACFS/R_Packages", .libPaths()))

#Clear everything that may be in the environment to start with

rm(list=ls())

# Load libraries

library(tidyverse)
library(data.table)
#library(meta)
library(devtools)
#library(pacman)
library(TwoSampleMR)
#library(MRInstruments)
library(ieugwasr)
#library(phenoscanner)
library(LDlinkR)
#library(MRPRESSO)
#library(extrafont)
#library(genetics.binaRies)

# Set working directory

setwd("/home/rmgpstc/Scratch/UKB")


# Run PLINK clumping from R

## First we will read in the exposure data from the relevant summary statistics file
# The below command is part of the TwoSampleMR package 
# From SSGAC Data Portal

exp_dat <- read_exposure_data(
  "EA4_additive_excl_23andMe.txt.gz",
  clump = FALSE,
  sep = "\t",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF_HRC",
  effect_allele_col = "Effect_allele",
  other_allele_col = "Other_allele",
  pval_col = "P",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "Chr",
  pos_col = "BP"
)

## Then we will clump instruments so that they are independent

# Clumping requires columns named rsid and pval so we rename these here

exp_dat <- exp_dat %>% 
  rename(
    rsid = SNP,
    pval = pval.exposure,
    )

# Then we use the ld_clump function in TwoSampleMR to clump.
# This step requires a token if the OpenGWAS is to be accessed.
# Instructions on how to get this can be found at 
# https://api.opengwas.io/api/#authentication or https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication
# Once the token is created, it can be saved in a .Renviron file in your home directory which can be found with command .libPaths()

exp_dat_clumped <- ld_clump(
  dat = exp_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  pop = "EUR",
 # #opengwas_jwt = get_opengwas_jwt(),  This is used if not using a .Renviron file
  bfile = "/home/rmgpstc/Scratch/UKB/EUR", #This is used if using a saved reference dataset on your computer rather than API
  plink_bin = "/home/rmgpstc/ACFS/Software/plink"
)

# We now rename these columns back to what they were for next stage

exp_dat_clumped <- exp_dat_clumped %>% 
  rename(
    SNP = rsid,
    pval.exposure = pval
  )

# NOw we pring the number of IVs for exposure

print(paste0("Number of IVs: ", as.character(length(exp_dat_clumped$SNP))))


# Now we can extract the outcome data using our clumped SNPs from the exposure

out_dat <- read_outcome_data(
  "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz",
  #  snps = exp_dat_clumped$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq_Tested_Allele_in_HRS",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  pval_col = "P",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "POS"
)


# Identify & print the exposure instruments missing from outcome GWAS

missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% out_dat$SNP)]
print(paste0("Number of IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
print("List of IVs missing from outcome GWAS:")
for (i in 1:length(missing_IVs)) {
  print(paste0(missing_IVs[i]))
}

# Now we replace the missing instruments from the outcome GWAS with proxies
# We do this using the package LDLinkR which automatically checks and replaces each rsID if a proxy is present
# This can take a bit of time

if(length(missing_IVs) == 0) {
  
  print("All exposure IVs found in outcome GWAS.")
  
} else {
  
  print("Some exposure IVs missing from outcome GWAS.")
  out_full <- fread("wmh_outcome_se.txt")
  
  for (i in 1:length(missing_IVs)) {
    
    proxies <- LDproxy(snp = missing_IVs[i], pop = "EUR", r2d = "r2", token = "79d1203b52c0", file = FALSE)
    proxies <- proxies[proxies$R2 > 0.8, ]
    proxy_present = FALSE
    
    if(length(proxies$RS_Number) == 0){
      
      print(paste0("No proxy SNP available for ", missing_IVs[i]))
      
    } else {
      
      for (j in 1:length(proxies$RS_Number)) {
        
        proxy_present <- proxies$RS_Number[j] %in% out_full$SNP
        
        if (proxy_present) {
          proxy_SNP = proxies$RS_Number[j]
          proxy_SNP_allele_1 = str_sub(proxies$Alleles[j], 2, 2)
          proxy_SNP_allele_2 = str_sub(proxies$Alleles[j], 4, 4)
          original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
          original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
          break
        }
      }
    }
    
    if(proxy_present == TRUE) {
      print(paste0("Proxy SNP found. ", missing_IVs[i], " replaced with ", proxy_SNP))
      proxy_row <- out_dat[1, ]
      proxy_row$SNP = missing_IVs[i]
      proxy_row$beta.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "Effect_beta"])
      proxy_row$se.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "se"])
      if (out_full[out_full$SNPID == proxy_SNP, "Allele1"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome = original_SNP_allele_1
      if (out_full[out_full$SNPID == proxy_SNP, "Allele1"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome = original_SNP_allele_2
      if (out_full[out_full$SNPID == proxy_SNP, "Allele2"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome = original_SNP_allele_1
      if (out_full[out_full$SNPID == proxy_SNP, "Allele2"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome = original_SNP_allele_2
      proxy_row$pval.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "pvalue"])
     # Blanking these next rows out as they're not in my dataset   
     # proxy_row$samplesize.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "N"])
     # if("N_case" %in% colnames(out_full)) proxy_row$ncase.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "N_case"])
     # if("N_control" %in% colnames(out_full))proxy_row$ncontrol.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "N_control"])
      proxy_row$chr.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "chr.exposure"])
      proxy_row$pos.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "pos.exposure"])
     # if("AF1" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNPID == proxy_SNP, "AF1"])
      out_dat <- rbind(out_dat, proxy_row)
    }
    
    if(proxy_present == FALSE) {
      print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
    }
  }
  
}

# Now that we have our exposure and outcome SNPs in the right format, we can harmonise

dat <- harmonise_data(
  exposure_dat = exp_dat_clumped, 
  outcome_dat = out_dat, 
  action = 1
)

saveRDS(dat, file = "datBMI.RData")

# And now we run our MR

datBMI <- readRDS("datBMI.RData")
res <- mr(datBMI, method_list = "mr_ivw")
