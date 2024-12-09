#!/bin/bash

# =============================================================================
# Title:        Running cutadapt for tRNA and rRNA removal
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

for i in ./cutadapt/SRR1*_trimmed.fastq.gz; do
    sample=${i##*/}
    sample=${sample%%_*}
    sbatch -o log/xrna_${sample}.out run_xrna.sh $i
done

