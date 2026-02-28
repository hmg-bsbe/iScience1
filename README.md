# iScience1 — RNA-seq Analysis Scripts

This repository contains the scripts used to perform standard differential gene expression (DGE) and transcriptome-wide analyses described in our manuscript - Gene expression signatures reveal accelerated brain ageing process in mouse models of neurodegenerative disorders.

No custom algorithms or novel computational methods were developed for this study. All analyses were performed using established and widely accepted RNA-seq workflows implemented in R. The scripts provided here represent the exact execution framework used to generate the results and figures reported in the manuscript.
---

## Scope of the Repository

The scripts implement standard analytical steps commonly used in RNA-seq studies, including:

* Import of gene-level abundance estimates
* Count normalization (e.g., TMM-based approaches)
* Transcriptome-wide exploratory analyses
* Correlation analysis of log₂ fold-change values across contrasts
* Generation of whole transcriptome heatmap, PCA and scatterplot figures

All statistical procedures rely on established Bioconductor and CRAN packages. No proprietary or unpublished pipelines were used.

---
## Repository Structure

**Scripts/**
Contains R scripts used for normalization, statistical testing, and figure generation. The scripts are modular and internally commented.

**Session_info/**
Contains R `sessionInfo()` outputs documenting software versions used during analysis to ensure reproducibility.
---

## Reproducibility Statement

The analyses were conducted using standard RNA-seq statistical frameworks. The repository provides the executable scripts used in the manuscript to enable full transparency and reproducibility of the reported results.

