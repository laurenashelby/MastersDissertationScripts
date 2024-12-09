# Probing Ribosome-mRNA Interactions: Start-Stop Elements vs. Canonical Sites

This repository contains the scripts, workflows, and visualizations developed as part of my Master's dissertation, *"Probing for differences in ribosome-mRNA interactions between start-stop elements and canonical start and stop positions."*

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation and Dependencies](#installation-and-dependencies)
- [Usage](#usage)
- [Data Availability](#data-availability)
- [Repository Structure](#repository-structure)
- [Acknowledgments](#acknowledgments)

## Overview

The project investigates start-stop (STST) elements in the 5' untranslated region (UTR) of mRNA, focusing on their interaction with ribosomes and comparing these interactions to canonical start and stop codons. Key methods include ribosome footprinting analysis of 40S and 80S subunits and advanced data visualization techniques.

## Features
- **Data Preprocessing:** Scripts for quality control, adapter trimming, alignment, and index building.
- **Visualization:** Heatmaps and density plots for ribosome footprint analysis.
- **Statistical Analysis:** Normalized counts and error-bar calculations for ribosomal positioning.

## Installation and Dependencies

### Dependencies
The following tools and packages are required:
- **Languages**: Python3, R
- **R Packages**: `tidyverse`, `FNN`, `DMwR2`
- **Bioinformatics Tools**: `FASTQC`, `Cutadapt`, `Bowtie`, `Bowtie2`, `STAR`, `Picard`, `SAMtools`
- **Database References**: `GENCODE43`, `UCSC Table Browser`

### Installation Steps
1. Clone this repository:
    ```bash
    git clone https://github.com/your-username/your-repo.git
    ```
2. Install the dependencies:
    - **Python**: Install packages via `pip install -r requirements.txt`.
    - **R**: Use `install.packages()` for listed libraries.
3. Set up bioinformatics tools as described in their documentation.

## Usage

### Data Preprocessing
Run scripts in the `Pre-Processing` directory sequentially to process raw ribosome footprint data:
1. **Download and Quality Control**: Use `sradownload.sh` to download the raw data from NCBI GEO and run `fastqc.sh` to perform quality control checks.
2. **Adapter Trimming**: Use `for_cutadapt.sh` to run `run_cutadapt_bioscientific.sh` for Bio-Scientific NEXTflex Small RNA-Seq Kit v3 libraries, and `run_cutadapt_takara.sh` for Takara/Clontech SMARTer smRNA-Seq libraries.
3. **tRNA and rRNA Removal**: Use `for_xrna.sh` for `run_xrna.sh`, with indices generated from `build-xrna-index.sbatch`.
4. **Alignment**: Build the genome index using `build-star-index.sbatch` with the human genome (`hg38.fa`) and gene annotations (`gencode.v43.chr_patch_hapl_scaff.annotation.gtf`), then use `for_star.sh` for `run_star.sh`.
5. **Parse and Compile Data**: Use `run_parse_sam.py`, then `run_parse_all_unique.sh` and `parse_processed_txt_strand.R` to generate TSV files for analysis.

### Visualization
Create heatmaps and density plots using R scripts:
- **Compile Data**: Use `1-transcript-count.Rmd` and `1.5-data-preparation.Rmd`.
- **Generate Visualizations**: Run `DissertationPlots.R` to access `plots_script_final.R` for heatmaps and density plots.

## Data Availability

Raw data utilized in this study are publicly available from the NCBI GEO database. Key datasets include:

- **GEO Accession Number**: [GSE140084](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140084)
- **80S Ribosome Reads**: `SRR10416857`
- **40S Ribosome Reads**: `SRR11553468`

## Repository Structure

```plaintext
Pre-Processing/
├── Index/
│   ├── build-star-index.sbatch
│   ├── build-xrna-index.sbatch
│   ├── gencode.v43.chr_patch_hapl_scaff.annotation.gtf
│   ├── hg38.fa
│   ├── human_rRNA.fa
│   └── human_tRNA.fa
├── SRA/
│   ├── fastqc.sh
│   ├── sradownload.sh
├── for_cutadapt.sh
├── run_cutadapt_bioscientific.sh
├── run_cutadapt_takara.sh
├── for_star.sh
├── run_star.sh
├── run_parse_sam.py
├── for_xrna.sh
├── run_xrna.sh
├── run_parse_all_unique.sh
├── parse_processed_txt_strand.R

Plots/
├── plots_script_final.R
├── 1-transcript-count.Rmd
├── 1.5-data-preparation.Rmd
└── DissertationPlots.Rmd
```

## Acknowledgments

This work was supported by mentorship and resources from the Vogel Lab at NYU. Special thanks to:

- **Dr. Christine Vogel**: Mentorship and guidance throughout the project.
- **Ziyue Cheng**: For providing foundational code and workflows.
- **Bohlen et al., 2020**: For their ribosome footprinting dataset that formed the basis of this analysis.



