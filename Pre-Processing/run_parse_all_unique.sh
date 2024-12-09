#!/bin/bash

# =============================================================================
# Title:        Conda Container for parse_processed_txt_strand.R
# Author:       Lauren Shelby, Ziyue Cheng
# Date:         December 9, 2024
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=200GB
#SBATCH --job-name=readInfo
##SBATCH --mail-type=END
#SBATCH --output=log/read_%j.out

######################################################
# INPUT: star/*_unique.txt.gz
# OUTPUT: parse-sam-double/*_with_read_info.tsv.gz
#     annotation with read info: length x position on transcript
#     e.g. 22x135,35x512
#
# Requirement: a singularity overlay with R installed by conda
#
# R packages: tidyverse, multidplyr
######################################################

singularity exec \
    -o /scratch/las10014/container/mapping.img:ro \
    /scratch/work/public/singularity/ubuntu-22.04.1.sif \
    /bin/bash -c "source /ext3/env.sh; conda activate base; Rscript ./parse_processed_txt_strand.R"

