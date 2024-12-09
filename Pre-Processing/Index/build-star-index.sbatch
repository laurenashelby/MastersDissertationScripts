#!/bin/bash

# =============================================================================
# Title:        Building STAR aligner index
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=STAR-build
#SBATCH --output=index_%j.out

module purge
module load star/intel/2.7.6a

STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir hg38_star \
    --genomeFastaFiles hg38.fa \
    --sjdbGTFfile gencode.v43.chr_patch_hapl_scaff.annotation.gtf \
    --sjdbOverhang 49