#!/bin/bash

# =============================================================================
# Title:        Run Quality Control on FASTQ files with FASTQC
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=fastqc
##SBATCH --mail-type=END
#SBATCH --output=fastqc_%j.out

module purge
module load fastqc/0.11.9

mkdir -p qc

for file in *.fastq.gz
do
    date
    fastqc -o qc $file
done

