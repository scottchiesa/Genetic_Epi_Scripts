#################################################################################################
###########################Performing a LiftOver Between Genome Builds###########################
####################################Author: Scott Chiesa 2025####################################
#################################################################################################

## Sometimes it is necessary to perform a liftover from a newer GRC Chromosome Build (GRCh38) to an older one (GRCh37) or vice versa.
## This ensures your coordinates align when comparing different datasets and/or reference panels and that the variants can therefore be correctly combined.

## The below example takes summary data from https://pubmed.ncbi.nlm.nih.gov/34012112/ and uses the UCSC Liftover tool to convert it from GRCh38 format back to GRCh37 (also known as hg19)


#==========Change into working directory==========#

cd path/to/your/directory


#=========Load necessary modules from Myriad==========#

module unload -f compilers mpi gcc-libs
module load r/recommended
module load python/miniconda3/24.3.0-0
source $UCL_CONDA_PATH/etc/profile.d/conda.sh


#==========Download necessary tools to perform liftover==========#

# First download Liftover tool from link below

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

# Then download the required chain file. This is the file that allows liftover to map the new and old co-ordinates.
# If doing the reverse liftover from GRCh37 to GRCh38 you'll need a different chain

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
chmod +x hg38ToHg19.over.chain.gz


#==========Download the summary stats file you want to liftover==========#

wget -O T1D.txt "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/GCST90014023_buildGRCh38.tsv"


#==========Prepare the file so it is in the correct format==========#

# Open R

R

# Read in file and then view first couple of rows to see how data are laid out

t1d_data <- fread("T1D.txt", header=T) 
head(t1d_data)

# Rename columns to standardised names based on order in file

names(t1d_data) <- c("RSID","P","CHR","BP","A1","A2","EAF","BETA","SE","N")

# Re-order columns so that CHR, BP, and RSID are in columns 1,2, and 3 respectively (rest can be any order).

t1d_data <- t1d_data[, c("CHR","BP","RSID","A1","A2","BETA","SE","EAF","P")]

# Write out to a txt file

fwrite(t1d_data, "Y:/UKB/Cristian Analysis/t1d_base_data.txt", sep = " ")


#==========Create virtual environment to run liftover==========#

## The UCSC Liftover tool requires libraries that Myriad doesn't have
## A conda environment with the liftover tool inside is therefore created so we run it from there instead

conda create -n liftover_env -c bioconda ucsc-liftover
conda activate liftover_env


#=========Prepare bash script that will be used to run command==========#

## The below code is a slightly modified version of that provided by CLS at https://github.com/CLS-Data/CLS_PGI_repository/blob/master/scripts/run_liftover.sh, so full credit to them for writing it.
## Save everything between the lines below as a bash file called run_analysis.sh

#====================================================================================================================================================================================================================================#

#!/bin/bash
# Convert GWAS summary stats from GRCh38 → GRCh37

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 [Input File] [Output File]"
    exit 1
fi

# If you are submitting this as a job to Myriad, the two lines below may need to be changed to be the names of your input and output files (with these then removed as arguments after the sh run_analysis line at the end)

INPUT_FILE="$1"
OUTPUT_FILE="$2"

LIFTOVER="liftOver"
CHAIN_FILE="hg38ToHg19.over.chain.gz"

if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Input file not found: ${INPUT_FILE}"
    exit 1
fi

BED_FILE="${INPUT_FILE%.txt}.bed"
UNMAPPED_FILE="${INPUT_FILE%.txt}_unmapped.txt"

echo "Creating BED file from ${INPUT_FILE}"
# adjust column numbers here to match your file structure! (this modified slightly from CLS code as nothing would map unless I set the end point (i.e. column 3) to be +1 the actual BP position due to coordinate system)
awk 'BEGIN {OFS="\t"} NR > 1 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {print "chr"$1, $2, $2+1, $3}' "${INPUT_FILE}" > "${BED_FILE}"


echo "Running LiftOver..."
"${LIFTOVER}" "${BED_FILE}" "${CHAIN_FILE}" "${OUTPUT_FILE}" "${UNMAPPED_FILE}"

if [[ -f "${OUTPUT_FILE}" ]]; then
    echo "✅ LiftOver completed: ${OUTPUT_FILE}"
else
    echo "❌ LiftOver failed"
fi

if [[ -s "${UNMAPPED_FILE}" ]]; then
    echo "Unmapped records saved: ${UNMAPPED_FILE}"
else
    echo "No unmapped variants."
fi

exit 0

#====================================================================================================================================================================================================================================#

## Now run the file you just saved to perform the liftover (the full command is 'sh [name of bash script] [name of input file] [name of output file]')

sh run_liftover.sh t1d_base_data.txt t1d_lifted_data.bed
