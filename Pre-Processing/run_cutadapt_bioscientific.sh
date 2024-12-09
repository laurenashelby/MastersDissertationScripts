#!/bin/bash

# =============================================================================
# Title:        Running cutadapt for Bioscientific specifications
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=cutadapt
#SBATCH --mail-type=END
#SBATCH --mail-user=las10014@nyu.edu
#SBATCH --output=log/cutadapt_%j.out


###############################
# Cut adapter of the sequence
#
# Toolkit of this paper: Bio-Scientific NEXTflex Small RNA-Seq Kit v3
# http://www.uusmb.unam.mx/assets/BibARN1.pdf
#
# -q quality > 10
# -a adapter sequence
# -m minimum length 10
# -e error rate 0.1
###############################

module purge
module load cutadapt/3.1

mkdir -p cutadapt

file=${1##*/}
file_name=${file%%.*}


cutadapt \
    -j $SLURM_CPUS_PER_TASK \
    -m 10 \
    -a TGGAATTCTCGGGTGCCAAGG \
    -q 10 \
    -e 0.1 \
    --discard-untrimmed \
    --discard-casava \
    -o cutadapt/${file_name}_trimmed.fastq.gz \
    $1
