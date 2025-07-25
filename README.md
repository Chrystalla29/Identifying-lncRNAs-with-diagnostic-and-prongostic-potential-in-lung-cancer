## LncRNA Biomarker Analysis in Lung Cancer (LUAD & LUSC)

This repository contains a comprehensive pipeline for identifying and analyzing long non-coding RNAs (lncRNAs) as diagnostic and prognostic biomarkers in two lung cancer subtypes: 
**LUAD (Lung Adenocarcinoma)** and **LUSC (Lung Squamous Cell Carcinoma)** using RNA-seq and clinical data from **The Cancer Genome Atlas (TCGA)**.

---

## OVERVIEW

The project includes:
- Downloading, processing and cleaning TCGA data 
- Performing differential expression analysis (DESeq2) on lncRNAs
- Data visualization
- Correlating expression levels with clino-demographic factors
- Kaplan-Meier survival analysis

---

## PROJECT CONTENT
- `data/`: Contains all CSV files for expression counts, differential expression, and survival summaries.
- `plots/`: Includes all visual outputs (e.g., MA plots, heatmaps, survival curves).
- `research_project_code.R`: The main R script for the full bioinformatics pipeline.
- `README.md`: You're here! The project overview and instructions.


## OBJECTIVE

To identify and evaluate differentially expressed lncRNAs in LUAD and LUSC as potential diagnostic and prognostic biomarkers and explore their association with clinical outcomes.

## METHODS OVERVIEW

- Data Source: TCGA RNA-Seq & clinical metadata (LUAD & LUSC)
- Tools: `TCGAbiolinks`, `DESeq2`, `pheatmap`, `ggplot2`, `survival`, `survminer`
- Pipeline:
  - Downloading RNA-seq and clinical metadata from TCGA
  - Extracting and filtering lncRNAs from raw gene expression data
  - Performing differential expression analysis (DESeq2)
  - Correlating expression levels with clinical and demographic variables
  - Visualizing results with MA plots, heatmaps, boxplots, etc
  - Conducting Kaplan-Meier survival analysis on top DE lncRNAs
  - Conducting Kaplan-Meier stratified survival analysis on top DE lncRNAs

## KEY OUTPUTS

- lncRNAs differentially expressed between tumor and normal samples
- lncRNAs with diagnostic potential in early stages of lung cancer
- lncRNAs with prognostic potential
- A "signature of lncRNAs", a group of lncRNAs that collectively characterize lung cancer.
  This signature is useful for patient stratification, dividing patients into groups like high-risk or low-risk based on their lncRNA patterns. 

## REQUIRMENTS

Install R libraries before running:

```r
BiocManager::install(c("TCGAbiolinks", "DESeq2", "pheatmap", "survival", "survminer", "ggplot2", "tibble", "dplyr", "RColorBrewer", "VennDiagram"))

## FUTURE WORK
-Further examine which lncRNAs can serve as diagnostic or prognostic tools or both
-Correlate DE lncRNAs with demographic and clinical features
