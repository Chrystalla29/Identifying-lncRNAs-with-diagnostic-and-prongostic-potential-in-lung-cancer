## LncRNA Biomarker Analysis in Lung Cancer (LUAD & LUSC)

This repository contains a comprehensive pipeline for identifying and analyzing long non-coding RNAs (lncRNAs) as diagnostic and prognostic biomarkers in two lung cancer subtypes: 
**LUAD (Lung Adenocarcinoma)** and **LUSC (Lung Squamous Cell Carcinoma)** using RNA-seq and clinical data from **The Cancer Genome Atlas (TCGA)**.

---

## OVERVIEW

The project includes:

- Downloading RNA-seq and clinical metadata from TCGA  
- Extracting and filtering lncRNAs from raw gene expression data  
- Performing differential expression analysis (DESeq2)  
- Visualizing results with MA plots, heatmaps, and boxplots  
- Conducting Kaplan-Meier survival analysis on top DE lncRNAs  

---

## PROJECT CONTENT
- `data/`: Contains all CSV files for expression counts, differential expression, and survival summaries.
- `plots/`: Includes all visual outputs (e.g., MA plots, heatmaps, survival curves).
- `research_project_code.R`: The main R script for the full bioinformatics pipeline.
- `README.md`: You're here! The project overview and instructions.


## OBJECTIVE

To identify and evaluate differentially expressed lncRNAs in LUAD and LUSC as potential diagnostic and prognostic biomarkers, and explore their association with clinical outcomes.

## METHODS OVERVIEW

- Data Source: TCGA RNA-Seq & clinical metadata (LUAD & LUSC)
- Tools: `TCGAbiolinks`, `DESeq2`, `pheatmap`, `ggplot2`, `survival`, `survminer`
- Pipeline:
  - Download and preprocess raw RNA-Seq data
  - Filter for lncRNAs using Ensembl annotations
  - Normalize counts and perform differential expression analysis (DEA)
  - Visualize DEGs via MA plots, heatmaps, and boxplots
  - Conduct survival analysis on top DEGs and significant lncRNAs

## KEY OUTPUTS

- MA plots: Differential expression landscape
- Boxplots: Top 5 DEGs per cancer type
- Heatmaps: Top 50 DE lncRNAs across samples
- Survival curves: Prognostic potential of top lncRNAs

## REQUIRMENTS

Install R libraries before running:

```r
BiocManager::install(c("TCGAbiolinks", "DESeq2", "pheatmap", "survival", "survminer", "ggplot2", "tibble", "dplyr", "RColorBrewer", "VennDiagram"))

## FUTURE WORK
-Further examine which lncRNAs can serve as diagnostic or prognostic tools or both
-Correlate DE lncRNAs with demographic and clinical features
