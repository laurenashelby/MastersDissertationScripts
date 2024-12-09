#!/bin/bash

# =============================================================================
# Title:        Running correct files for Takara and Bioscientific cutadapt
# Author:       Lauren A. Shelby
# Date:         December 9, 2024
# =============================================================================


for i in SRA/SRR10*.fastq.gz; do
    sbatch ./run_cutadapt_bioscientific.sh $i
done

for i in SRA/SRR11*.fastq.gz; do
    sbatch ./run_cutadapt_takara.sh $i
done
