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
- `*.R` : Scripts for data preprocessing, differential expression analysis, survival analysis, and correlation between lncRNAs and mRNAs.
- `sessionInfo.txt` : R environment details to reproduce the analysis.
- `README.md` : Instructions and project overview.
-  `plots/` : Figures and visual outputs (MA plots, PCA, boxplots, survival curves, heatmaps).


## OBJECTIVE

To identify and evaluate differentially expressed lncRNAs in LUAD and LUSC as potential diagnostic and prognostic biomarkers and explore their association with clinical outcomes.

## DATA
Raw RNA-seq and clinical data were downloaded from [TCGA](https://portal.gdc.cancer.gov/).  
Data files are **not included** in this repository due to size limitations.  

## ENVIRONMENT
This project was developed in R version 4.4.2 (Windows 11).  
All package versions are listed in [sessionInfo.txt](sessionInfo.txt).

To reproduce the environment, install the required packages:
```r
# Example installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment",
                       "GenomicRanges",
                       "TCGAbiolinks",
                       "DESeq2",
                       "survminer",
                       "corrplot"))
```
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

