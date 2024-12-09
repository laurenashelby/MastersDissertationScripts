#!/bin/bash

# =============================================================================
# Title:        Running STAR aligner
# Author:       Lauren A. Shelby
# Date:         December 9, 2024
# =============================================================================

######################################################
# get mapped results
#
# INPUT: xrna/*_trimmed_xrna.fastq.gz
# OUTPUT: 
#   - hisat2/*.bam
#   - hisat2/*_summary.txt
#   - hisat2/*_unique.txt.gz
######################################################

for i in ./xrna/SRR1*_trimmed_xrna.fastq
do
    sbatch ./run_star.sh $i
done
