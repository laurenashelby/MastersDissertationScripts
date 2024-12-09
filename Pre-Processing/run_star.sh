#!/bin/bash

# =============================================================================
# Title:        Running STAR aligner, converting SAM to BAM, then SAM to txt
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=star
##SBATCH --mail-type=END
#SBATCH --output=log/star_%j.out

###############################
# Alignment with STAR
#
# OUTPUT: 
#   - star/*.Aligned.bam
#   - star/*_unique.txt.gz
###############################

echo [`date +"%m-%d-%y %T"`] START with job ID: $SLURM_JOBID
START_TIME=`date +%s`

mkdir -p star

file=${1##*/}
file_name=${file%.fastq*}
file_name=${file_name%%_*}

module purge
module load star/intel/2.7.11a

if [[ $1 =~ gz$ ]]; then
    READCOMMAND="--readFilesCommand zcat" # if it is a gzip file, read with zcat
fi

STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir ./index/hg38_star/ \
    --alignEndsType EndToEnd \
    --readFilesIn $1 \
    $READCOMMAND \
    --outFileNamePrefix ./star/${file_name}. \
    --outSAMtype BAM SortedByCoordinate

module load picard/2.27.5

java -jar $PICARD_JAR MarkDuplicates \
    I=./star/${file_name}.Aligned.sortedByCoord.out.bam \
    O=./star/${file_name}.Aligned.sorted.marked.noDups.bam \
    M=./star/${file_name}_marked_dup_metrics.txt \
    REMOVE_DUPLICATES=true

# get unique reads
echo [`date +"%m-%d-%y %T"`] Unique Reads

module load samtools/intel/1.14
module load python/intel/3.8.6

#    awk -v OFS='\t' '{print $1,$2,$3,$4,$6,length($10)}' | \
samtools view \
    -@ $SLURM_CPUS_PER_TASK \
    -d NH:1 \
    ./star/${file_name}.Aligned.sorted.marked.noDups.bam | \
    python ./run_parse_sam.py | \
        pigz -p $SLURM_CPUS_PER_TASK \
        > star/${file_name}_unique.txt.gz

END_TIME=`date +%s`
RUNNING_TIME=$((END_TIME - START_TIME))
echo [`date +"%m-%d-%y %T"`] END with job ID: $SLURM_JOBID, Elapsed: $(($RUNNING_TIME / 3600))hrs $((($RUNNING_TIME / 60) % 60))min $(($RUNNING_TIME % 60))sec

