#!/bin/bash

# =============================================================================
# Title:        Building Bowtie and Bowtie2 indices for rRNA and tRNA removal
# Author:       Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

module purge 
module load bowtie2/2.4.4

# bowtie rRNA
module purge
module load bowtie/1.3.1
bowtie-build human_rRNA.fa human_rRNA

# bowtie2 tRNA
module purge
module load bowtie2/2.4.4
bowtie2-build human_tRNA.fa human_tRNA
