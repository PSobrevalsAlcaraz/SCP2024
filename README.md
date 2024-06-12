# Single Cell Proteomics (SCP) Analysis Workflow

This repository contains an R script that executes the complete workflow for single cell proteomics analysis. The script is designed to be user-friendly, with an intuitive interface, and comes with a PDF document that provides a detailed explanation of each step.
Adapted from the methodology by Gatto and Vanderaa (2021 & 2023).

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow Description](#workflow-description)
- [Reference](#reference)
- [Contributing](#contributing)

## Features

- Comprehensive workflow for single cell proteomics analysis
- User-friendly interface
- Detailed documentation for each step
- Optimised for Proteome discoverer output data
- Adapted from the methodology by Gatto and Vanderaa (2021 & 2023)

## Installation

To use this workflow, you need to have R installed on your system. You can download R from [CRAN](https://cran.r-project.org/), this script has been made with R (version 4.1.2).
You should also have the following dependencies installed:
- devtools  - scp  - limma  - scuttle  - matrixStats  - GGally  - randomcoloR  - expss  - svDialogs  - dplyr  - scater  - ggbreak  - viridis  - seqinr  - ggpubr  - ggplot2  - RColorBrewer  - scp  - magrittr  - tidyr  - proDA  - tcltk  - anchors  - ggrepel  - plotly  - DEP  - reshape2  - ggthemes

If not, the script will install them automatically at the beginning.

## Usage

1. Open the R script `SCP_2024.R` in your R environment.
2. Follow the instructions provided in the script to load your data and run the analysis.
3. The input files provided for the analysis should be the output files from Proteome Discoverer (currenty based on PD 3.0):
  A) PSM.txt containing the following columns: File.ID, Master.Protein.Accessions, Annotated.Sequence, Percolator.PEP, Abundance."TMT". Only normal abundances should be selected; no other normalized abundance columns should be added.
  B) InputFiles.txt
5. Refer to the PDF document `SCP_workflow_2024.pdf` for a detailed explanation of each step.

## Workflow Description

The workflow is divided into several steps, listed below:

1. **PSMs (Peptide-Spectrum Matches)**: Files filtering according to # PSMs.
2. **Median SCR (Single Cell Ratio)**: Compute the Median Single Cell Ratio against the booster to filter out unexpected PSMs.
3. **FDR (False Discovery Rate)**: Compute FDR as q-values from posterior error probabilities (PEPs)
4. **Peptides**: Group PSMs to peptides using the median of PSMs.
5. **Median RI (Relative Intensity)**: Compute Median Relative reporter ion Intensity to filter low-quality cells.
6. **Median CV (Coefficient of Variance)**: Determine the variability of protein intensities.
7. **Number of Unique Peptides**: Filter proteins by the number of unique peptides
8. **Proteins**: Group peptides into proteins using the median of peptides.
9. **QC Metrics (Quality Control)**: Calculate QC metrics to identify and remove low-quality cells.
10. **Missingness Percentage**: Compute the percentage of missingness for each protein, and decide threshold of imputation.
12. **Imputation**: Impute missing values using random sampling or KNN.
13. **Batch Correction**: Perform batch correction on chosen values Samples or TMT.
14. **Gene Names**: Translate protein accession numbers to gene names.
15. **PCA (Principal Component Analysis)**
16. **UMAP (Uniform Manifold Approximation and Projection)**
    
Please refer to the included PDF document `SCP_workflow_2024.pdf` for a detailed description of each step.

## Reference

This workflow is based on the methodology described by Laurent Gatto and Christophe Vanderaa in their publications:

- Vanderaa, Christophe, and Laurent Gatto. 2021. “Replication of Single-Cell Proteomics Data Reveals Important Computational Challenges.” Expert Rev. Proteomics, October.
- Gatto, L., & Vanderaa, C. (2024). Single Cell Proteomics data processing and analysis (Version 1.14.0) Retrieved from [bioconductor](https://bioconductor.org/packages/release/bioc/vignettes/scp/inst/doc/scp.html)

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request.
