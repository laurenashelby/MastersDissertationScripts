#!/bin/bash

# =============================================================================
# Title:        Download FASTQ files with SRA tools
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=fq-dump
##SBATCH --mail-type=END
#SBATCH --output=fa-dump_%j.out

#################################
# download all the homo sapien  


module purge
module load sra-tools/3.1.0

SRR_ACC=("SRR10416857" "SRR11553468")

for i in "${SRR_ACC[@]}"
do
    echo $i
    fastq-dump --gzip $i
done
