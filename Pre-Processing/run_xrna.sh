#!/bin/bash

# =============================================================================
# Title:        Running Bowtie and Bowtie2 for tRNA and rRNA removal
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=bowtie
##SBATCH --mail-type=END
#SBATCH --output=log/xrna_%j.out

###############################
# Remove rRNA and tRNA using the Bowtie and Bowtie2
#
# with --un parameter, you can save the sequence with out rRNA or tRNA
#
# OUTPUT: xrna/*_xran.fastq.gz
###############################

echo [`date +"%m-%d-%y %T"`] START with job ID: $SLURM_JOBID
START_TIME=`date +%s`

file=${1##*/}
file_name=${file%%.*}

mkdir -p xrna

module purge
module load bowtie/1.3.1

echo [`date +"%m-%d-%y %T"`] Removing rRNA with Bowtie

bowtie \
    --threads $SLURM_CPUS_PER_TASK \
    --seedlen 18 \
    --un ./xrna/${file_name}_no_rRNA.fastq \
    -x ./index/human_rRNA \
    $1 > /dev/null

module purge
module load bowtie2/2.4.4

echo [`date +"%m-%d-%y %T"`] Removing tRNA with Bowtie2

bowtie2 \
    --threads $SLURM_CPUS_PER_TASK \
    -L 18 \
    --un ./xrna/${file_name}_xrna.fastq \
    -x ./index/human_tRNA \
    -U ./xrna/${file_name}_no_rRNA.fastq > /dev/null

rm ./xrna/${file_name}_no_rRNA.fastq

echo [`date +"%m-%d-%y %T"`] Compressing the fastq

pigz -p $SLURM_CPUS_PER_TASK ./xrna/${file_name}_xrna.fastq

END_TIME=`date +%s`
RUNNING_TIME=$((END_TIME - START_TIME))
echo [`date +"%m-%d-%y %T"`] END with job ID: $SLURM_JOBID, Elapsed: $(($RUNNING_TIME / 3600))hrs $((($RUNNING_TIME / 60) % 60))min $(($RUNNING_TIME % 60))sec

