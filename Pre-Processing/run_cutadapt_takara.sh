#!/bin/bash

# =============================================================================
# Title:        Running cutadapt for Takara specifications
# Author:       Lauren A. Shelby
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
#
# Toolkit of this paper: Bio-Scientific NEXTflex Small RNA-Seq Kit v3
# https://www.takarabio.com/documents/User%20Manual/SMARTer%20smRNA/SMARTer%20smRNA-Seq%20Kit%20for%20Illumina%20User%20Manual.pdf
#
#
# -q quality > 10 from 3' end
# -O overlay between the read and adapter
# -m minimum length 15
# -u trim first 3 nucleotides
# -e error rate 0.1
###############################

module purge
module load cutadapt/3.1

mkdir -p cutadapt

file=${1##*/}
file_name=${file%%.*}

cutadapt \
    -j $SLURM_CPUS_PER_TASK \
    -m 15 \
    -u 3 \
    -a AAAAAAAAAA \
    -q 10 \
    -e 0.1 \
    --discard-untrimmed \
    --discard-casava \
    -o cutadapt/${file_name}_trimmed.fastq.gz \
    $1

