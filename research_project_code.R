### RESEARCH PROJECT
#----------------------------------------------------------------------------------

#----------------DOWNLOADING DATA FOR LUAD AND LUSC-----------------------------
# Loading required libraries
library(TCGAbiolinks)  # Bioconductor package for querying and downloading TCGA data

# Get a list of available TCGA projects
gdcprojects <- getGDCprojects()

# Print summary information about the selected projects (LUAD and LUSC)
getProjectSummary("TCGA-LUAD")  # Lung Adenocarcinoma (LUAD)
getProjectSummary("TCGA-LUSC")  # Lung Squamous Cell Carcinoma (LUSC)


# Query for LUAD RNA-Seq data
query_LUAD <- GDCquery(
  project = "TCGA-LUAD",  # TCGA project ID
  data.category = "Transcriptome Profiling",  # Category for RNA sequencing data
  data.type = "Gene Expression Quantification",  # Retrieve counts per gene
  experimental.strategy = "RNA-Seq",  # Ensures only RNA-Seq data is retrieved
  workflow.type = "STAR - Counts"  # STAR alignment workflow for count data
)

# Query for LUSC RNA-Seq data
query_LUSC <- GDCquery(
  project = "TCGA-LUSC",  
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts"
)

# Download the RNA-Seq data for LUAD and LUSC
GDCdownload(query_LUAD)
GDCdownload(query_LUSC)

# Prepare the downloaded data for analysis 
# converts raw GDC format to summarized experiment object containing:
# gene expression matrix, sample metadata, gene metadata
LUAD_data <- GDCprepare(query_LUAD)
LUSC_data <- GDCprepare(query_LUSC)

#---------------------------DOWNLOAD lncRNA DATA-------------------------------
# Installing and loading annotation packages
#BiocManager::install("ensembldb")  # Package for working with Ensembl databases
library(ensembldb)
#BiocManager::install("GenomicRanges")  # Provides genomic range operations
library(GenomicRanges)
#BiocManager::install("GenomeInfoDb")  # Required for working with genome annotations
library(GenomeInfoDb)
#BiocManager::install("EnsDb.Hsapiens.v86")  # Human genome annotation database
library(EnsDb.Hsapiens.v86)

# Check all available genes from the Ensembl annotation database
all_genes <- genes(EnsDb.Hsapiens.v86)

# Check available gene biotypes to identify relevant categories for lncRNAs
unique(all_genes$gene_biotype)

# Define a list of biotypes corresponding to long non-coding RNAs (lncRNAs)
lncRNA_types <- c("lincRNA", "macro_lncRNA", "3prime_overlapping_ncRNA", 
                  "bidirectional_promoter_lncRNA", "processed_transcript", 
                  "sense_intronic", "sense_overlapping", "non_coding")


# Extract only the genes classified as lncRNAs
lncRNA_genes <- all_genes[all_genes$gene_biotype %in% lncRNA_types, ]

# Print a few lncRNA gene IDs to verify extraction
head(lncRNA_genes$gene_id)

# Ensure gene IDs in RNA-seq data match Ensembl format 
rownames(LUAD_data) <- sub("\\..*", "", rownames(LUAD_data))
rownames(LUSC_data) <- sub("\\..*", "", rownames(LUSC_data))

#MATCHING lncRNAs IN DATA
# Find lncRNA genes that are present in LUAD RNA-Seq data
lncRNA_ids_LUAD <- intersect(rownames(LUAD_data), lncRNA_genes$gene_id)

# Find lncRNA genes that are present in LUSC RNA-Seq data
lncRNA_ids_LUSC <- intersect(rownames(LUSC_data), lncRNA_genes$gene_id)

# Check the number of matched lncRNAs in both datasets
length(lncRNA_ids_LUAD) 
length(lncRNA_ids_LUSC)  

# Extract only lncRNA genes from LUAD and LUSC datasets
LUAD_lncRNA <- LUAD_data[lncRNA_ids_LUAD, ]
LUSC_lncRNA <- LUSC_data[lncRNA_ids_LUSC, ]

#SAVING FILTERED DATA
# Convert the extracted SummarizedExperiment objects to data frames for saving
LUAD_counts_df <- as.data.frame(assay(LUAD_lncRNA))
LUSC_counts_df <- as.data.frame(assay(LUSC_lncRNA))

# Save the extracted lncRNA expression data as CSV files
write.csv(LUAD_counts_df, "LUAD_lncRNA.csv", row.names = TRUE)
write.csv(LUSC_counts_df, "LUSC_lncRNA.csv", row.names = TRUE)

#-----------------------DOWNLOADING FOR NORMAL SAMPLES--------------------------
# Query LUAD including Normal Samples
query_LUAD_control <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"  # Ensuring we get normal samples
)

# Query LUSC including Normal Samples
query_LUSC_control <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

# Download normal sample data
GDCdownload(query_LUAD_control)
GDCdownload(query_LUSC_control)

# Prepare data for analysis
LUAD_control_data <- GDCprepare(query_LUAD_control)
LUSC_control_data <- GDCprepare(query_LUSC_control)

# Remove version numbers from gene IDs
rownames(LUAD_control_data) <- sub("\\..*", "", rownames(LUAD_control_data))
rownames(LUSC_control_data) <- sub("\\..*", "", rownames(LUSC_control_data))

# Find lncRNA genes present in normal samples
lncRNA_ids_LUAD_control <- intersect(rownames(LUAD_control_data), lncRNA_genes$gene_id)
lncRNA_ids_LUSC_control <- intersect(rownames(LUSC_control_data), lncRNA_genes$gene_id)

# Extract only lncRNA genes from normal samples
LUAD_lncRNA_control <- LUAD_control_data[lncRNA_ids_LUAD_control, ]
LUSC_lncRNA_control <- LUSC_control_data[lncRNA_ids_LUSC_control, ]

# Convert to data frames
LUAD_control_counts_df <- as.data.frame(assay(LUAD_lncRNA_control))
LUSC_control_counts_df <- as.data.frame(assay(LUSC_lncRNA_control))

# Save the extracted lncRNA expression data as CSV files
write.csv(LUAD_control_counts_df, "LUAD_control_lncRNA.csv", row.names = TRUE)
write.csv(LUSC_control_counts_df, "LUSC_control_lncRNA.csv", row.names = TRUE)


#----------------------Log transform normal counts -----------------------------
#Load Required Libraries 
library(DESeq2)
# Read raw counts
LUAD_control_counts <- read.csv("LUAD_control_lncRNA.csv", row.names = 1)
LUSC_control_counts <- read.csv("LUSC_control_lncRNA.csv", row.names = 1)

# Make dummy colData
LUAD_control_colData <- data.frame(condition = rep("Normal", ncol(LUAD_control_counts)))
rownames(LUAD_control_colData) <- colnames(LUAD_control_counts)

LUSC_control_colData <- data.frame(condition = rep("Normal", ncol(LUSC_control_counts)))
rownames(LUSC_control_colData) <- colnames(LUSC_control_counts)

# Create DESeq2 object and VST
dds_LUAD_control <- DESeqDataSetFromMatrix(countData = LUAD_control_counts, colData = LUAD_control_colData, design = ~ 1)
dds_LUSC_control <- DESeqDataSetFromMatrix(countData = LUSC_control_counts, colData = LUSC_control_colData, design = ~ 1)

# Run VST
LUAD_control_vst <- vst(dds_LUAD_control, blind = TRUE)
LUSC_control_vst <- vst(dds_LUSC_control, blind = TRUE)

# Save to CSV
write.csv(assay(LUAD_control_vst), "LUAD_control_lncRNA_log_trans_after_norm.csv")
write.csv(assay(LUSC_control_vst), "LUSC_control_lncRNA_log_trans_after_norm.csv")

#---------------------------------NORMALIZE COUNTS------------------------------
#Load Required Libraries 
library(DESeq2)
library(tibble)
library(readr)

# Loading data
# Read in tumor and control lncRNA counts for LUAD and LUSC
LUAD_counts <- read.csv("LUAD_lncRNA.csv", row.names = 1)
summary(LUAD_counts)
LUAD_controls <- read.csv("LUAD_control_lncRNA.csv", row.names = 1)

LUSC_counts <- read.csv("LUSC_lncRNA.csv", row.names = 1)
LUSC_controls <- read.csv("LUSC_control_lncRNA.csv", row.names = 1)

# Fix sample names (TCGA usually uses '-' not '.') so samples will match
fix_colnames <- function(df) {
  colnames(df) <- gsub("\\.", "-", colnames(df))
  return(df)
}

LUAD_counts <- fix_colnames(LUAD_counts)
LUAD_controls <- fix_colnames(LUAD_controls)
LUSC_counts <- fix_colnames(LUSC_counts)
LUSC_controls <- fix_colnames(LUSC_controls)

# Add suffix to each sample to distinguish the tumor and normal samples
make_unique_colnames <- function(df, label) {
  colnames(df) <- paste0(gsub("\\.", "-", colnames(df)), "_", label)
  return(df)
}

LUAD_counts <- make_unique_colnames(LUAD_counts, "tumor")
LUAD_controls <- make_unique_colnames(LUAD_controls, "normal")
LUSC_counts <- make_unique_colnames(LUSC_counts, "tumor")
LUSC_controls <- make_unique_colnames(LUSC_controls, "normal")


# Combine tumor and control counts for LUAD
LUAD_combined <- cbind(LUAD_counts, LUAD_controls)
LUAD_condition <- factor(c(rep("Tumor", ncol(LUAD_counts)), rep("Normal", ncol(LUAD_controls))))

# Combine tumor and control counts for LUSC
LUSC_combined <- cbind(LUSC_counts, LUSC_controls)
LUSC_condition <- factor(c(rep("Tumor", ncol(LUSC_counts)), rep("Normal", ncol(LUSC_controls))))

# Save for future use if needed
write.csv(LUAD_combined, "LUAD_combined.csv")
write.csv(LUSC_combined, "LUSC_combined.csv") 

#Get the number of samples
cat("LUAD Tumor samples:", ncol(LUAD_counts), "\n")
cat("LUAD Normal samples:", ncol(LUAD_controls), "\n")
cat("LUSC Tumor samples:", ncol(LUSC_counts), "\n")
cat("LUSC Normal samples:", ncol(LUSC_controls), "\n")


# Sample metadata for DESeq2 (first column = sample, second column=tumor or normal)
LUAD_colData <- data.frame(condition = LUAD_condition, row.names = colnames(LUAD_combined))
LUSC_colData <- data.frame(condition = LUSC_condition, row.names = colnames(LUSC_combined))

#Create DESeq2 Objects 
dds_LUAD <- DESeqDataSetFromMatrix(countData = LUAD_combined,
                                   colData = LUAD_colData,
                                   design = ~ condition)

dds_LUSC <- DESeqDataSetFromMatrix(countData = LUSC_combined,
                                   colData = LUSC_colData,
                                   design = ~ condition)

# Filter lowly expressed genes 
#Keeping geneswith at least 10 reads in 5 or more samples
keep_LUAD <- rowSums(counts(dds_LUAD) >= 10) >= 5
keep_LUSC <- rowSums(counts(dds_LUSC) >= 10) >= 5
dds_LUAD <- dds_LUAD[keep_LUAD, ]
dds_LUSC <- dds_LUSC[keep_LUSC, ]

# Run DESeq2 for normalization
dds_LUAD <- DESeq(dds_LUAD)
dds_LUSC <- DESeq(dds_LUSC)

# Save DESeq2 object after creating it
saveRDS(dds_LUAD, file = "dds_LUAD.rds")
saveRDS(dds_LUSC, file = "dds_LUSC.rds")

#Get Normalized Counts 
LUAD_norm_counts <- counts(dds_LUAD, normalized = TRUE)
LUSC_norm_counts <- counts(dds_LUSC, normalized = TRUE)


#Save normalized counts
write.csv(LUAD_norm_counts, "LUAD_lncRNA_normalized_counts_DESeq2.csv")
write.csv(LUSC_norm_counts, "LUSC_lncRNA_normalized_counts_DESeq2.csv")

# Log Transform for visualization 
LUAD_vst <- vst(dds_LUAD, blind = TRUE)
LUSC_vst <- vst(dds_LUSC, blind = TRUE)

#Save the log transformed data
write.csv(assay(LUAD_vst), "LUAD_lncRNA_log_trans_after_norm.csv")
write.csv(assay(LUSC_vst), "LUSC_lncRNA_log_trans_after_norm.csv")

#=====Before and after VST plot and normalization for LUAD and LUSC=============

# Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)

#Function for processing and plotting
plot_expression <- function(raw_file, norm_file, vst_file, cancer_type) {
  
  # Load data
  raw <- read.csv(raw_file, row.names = 1)
  norm <- read.csv(norm_file, row.names = 1)
  vst <- read.csv(vst_file, row.names = 1)
  
  # Keep only tumor samples
  tumor_samples <- grep("_tumor$", colnames(raw), value = TRUE)
  
  # Select 30 random tumor samples
  set.seed(310)
  selected_samples <- sample(tumor_samples, 30)
  
  # Subset data
  raw_subset <- raw[, selected_samples]
  vst_subset <- vst[, selected_samples]
  
  # Melt to long format
  raw_long <- melt(raw_subset, variable.name = "Sample", value.name = "Counts")
  vst_long <- melt(vst_subset, variable.name = "Sample", value.name = "Counts")
  
  # Label stage
  raw_long$Stage <- "Raw"
  vst_long$Stage <- "VST"
  
  # Combine and order
  combined <- rbind(raw_long, vst_long)
  combined$Stage <- factor(combined$Stage, levels = c("Raw", "VST"))
  
  # Remove extreme outliers
  combined <- combined %>% filter(Counts < quantile(Counts, 0.999))
  
  # Shorten sample labels and add group
  short_id <- sapply(strsplit(as.character(combined$Sample), "\\."), function(x) paste(x[1:3], collapse = "-"))
  group <- ifelse(grepl("tumor$", combined$Sample, ignore.case = TRUE), "Tumor", "Normal")
  combined$Sample <- paste0(short_id, " (", group, ")")
  
  # Plot
  ggplot(combined, aes(x = Sample, y = Counts, fill = Stage)) +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~Stage, nrow = 1, scales = "free_y") +
    theme_bw(base_size = 14) +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y  = element_text(size = 14),
      strip.text   = element_text(size = 16, face = "bold"),
      plot.title   = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    ) +
    labs(title = paste0(cancer_type, ": Raw vs VST Expression"),
         y = "Expression", x = "Sample") +
    scale_fill_manual(values = c("Raw" = "red", "VST" = "#F4A8A8"))
}

#Generate LUAD and LUSC plots

# LUAD
plot_expression(
  raw_file = "LUAD_combined.csv",
  norm_file = "LUAD_lncRNA_normalized_counts_DESeq2.csv",
  vst_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  cancer_type = "LUAD"
)

# LUSC
plot_expression(
  raw_file = "LUSC_combined.csv",
  norm_file = "LUSC_lncRNA_normalized_counts_DESeq2.csv",
  vst_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  cancer_type = "LUSC"
)





#----------------------------- Differential Expression -------------------------
#Perform DEA
# These results contain:
# - log2FoldChange: expression difference between Tumor and Normal
# - pvalue: raw p-value
# - padj: adjusted p-value (FDR)
# - baseMean: average expression across all samples
res_LUAD <- results(dds_LUAD, contrast = c("condition", "Tumor", "Normal"))
res_LUSC <- results(dds_LUSC, contrast = c("condition", "Tumor", "Normal"))

#Shrink log2 fold changes
library(apeglm)
res_LUAD_shrink <- lfcShrink(dds_LUAD,coef = "condition_Tumor_vs_Normal",type = "apeglm")  
res_LUSC_shrink <- lfcShrink(dds_LUSC,coef = "condition_Tumor_vs_Normal",type = "apeglm")  

#Filter results by significance
#Keep genes with: FDR-adjusted p-value (padj) < 0.01 and absolute log2 Fold Change > 1 
res_LUAD_filtered <- subset(res_LUAD_shrink, padj < 0.01 & abs(log2FoldChange) > 1)
res_LUSC_filtered <- subset(res_LUSC_shrink, padj < 0.01 & abs(log2FoldChange) > 1)

# Save full and filtered results
write.csv(as.data.frame(res_LUAD), "LUAD_DESeq2_all_results.csv")
write.csv(as.data.frame(res_LUSC), "LUSC_DESeq2_all_results.csv")
write.csv(as.data.frame(res_LUAD_filtered), "LUAD_DESeq2_all_results_shrunk.csv")
write.csv(as.data.frame(res_LUSC_filtered), "LUSC_DESeq2_all_results_shrunk.csv")
#---------------------------------MA plot---------------------------------------
#Load library
library(tidyverse)

# Load DESeq2 results 
res_LUAD <- read.csv("LUAD_DESeq2_all_results.csv", row.names = 1)
res_LUSC <- read.csv("LUSC_DESeq2_all_results.csv", row.names = 1)

####For LUAD
p <- ggplot(res_LUAD, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "#F4A8A8", "TRUE" = "#A80000")) + # light red for non-sig, deep red for sig
  scale_x_log10() +
  labs(
    title = "MA Plot - LUAD ",
    x = "Mean Expression (log10)",
    y = "log2 Fold Change",
    color = "padj < 0.05"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),    # Title size
    axis.title.x = element_text(size = 24),                              # X-axis label size
    axis.title.y = element_text(size = 24),                              # Y-axis label size
    axis.text = element_text(size = 24),                                 # Axis tick label size
    legend.title = element_text(size = 24),                              # Legend title
    legend.text = element_text(size = 24),                               # Legend text
    legend.position = "bottom"
  )

#Save plot
ggsave("MA_plot_LUAD.png", plot = p, width = 8, height = 6, dpi = 300)


####For LUSC
p2 <- ggplot(res_LUSC, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "#F4A8A8", "TRUE" = "#A80000")) + 
  scale_x_log10() +
  labs(
    title = "MA Plot - LUSC ",
    x = "Mean Expression (log10)",
    y = "log2 Fold Change",
    color = "padj < 0.05"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),    # Title size
    axis.title.x = element_text(size = 24),                              # X-axis label size
    axis.title.y = element_text(size = 24),                              # Y-axis label size
    axis.text = element_text(size = 24),                                 # Axis tick label size
    legend.title = element_text(size = 24),                              # Legend title
    legend.text = element_text(size = 24),                               # Legend text
    legend.position = "bottom"
  )


#Save the plot
ggsave("MA_plot_LUSC.png", plot = p2, width = 8, height = 6, dpi = 300)



#=======MA plots highlighting also fold change=====
#Load required library
library(tidyverse) 

# Load DESeq2 results
# These CSV files should be generated by DESeq2 and contain columns like baseMean, log2FoldChange, and padj
res_LUAD <- read.csv("LUAD_DESeq2_all_results.csv", row.names = 1)
res_LUSC <- read.csv("LUSC_DESeq2_all_results.csv", row.names = 1)

#Function to classify gene significance and direction
classify_significance <- function(df) {
  df$significance <- "Not Sig"  # Default classification
  df$significance[df$padj < 0.05 & df$log2FoldChange > 0] <- "Up"    # Significantly upregulated
  df$significance[df$padj < 0.05 & df$log2FoldChange < 0] <- "Down"  # Significantly downregulated
  df$significance <- factor(df$significance, levels = c("Up", "Down", "Not Sig"))  # Set factor levels
  return(df)
}

# Apply classification to both datasets 
res_LUAD <- classify_significance(res_LUAD)
res_LUSC <- classify_significance(res_LUSC)

#Define color palette 
color_palette <- c(
  "Up" = "#A80000",         # Deep red for upregulated
  "Down" = "#F4A8A8",       # Light red for downregulated
  "Not Sig" = "grey"     # Grey for non-significant
)

#Create MA plot for LUAD
plot_LUAD <- ggplot(res_LUAD, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significance), alpha = 0.6) +  # Scatter plot with alpha for transparency
  scale_color_manual(values = color_palette) +          # Apply custom colors
  scale_x_log10() +                                     # Log-scale for x-axis
  labs(
    title = "MA Plot - LUAD",                           # Plot title
    x = "Mean Expression (log10)",                      # X-axis label
    y = "log2 Fold Change",                             # Y-axis label
    color = "Gene Regulation"                           # Legend title
  ) +
  theme_minimal(base_size = 14) +                       # Clean theme with readable font
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),  # Centered bold title
    legend.position = "bottom"                              # Legend at bottom
  )

# Create MA plot for LUSC
plot_LUSC <- ggplot(res_LUSC, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(values = color_palette) +
  scale_x_log10() +
  labs(
    title = "MA Plot - LUSC",
    x = "Mean Expression (log10)",
    y = "log2 Fold Change",
    color = "Gene Regulation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

#Save plots to PNG files 
ggsave("MA_Plot_LUAD2.png", plot_LUAD, width = 8, height = 6, dpi = 300)
ggsave("MA_Plot_LUSC2.png", plot_LUSC, width = 8, height = 6, dpi = 300)


print(plot_LUAD)
print(plot_LUSC)



#------------------------HEATMAP using expression values------------------------
# Load required libraries
library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Define Heatmap Function for Top 50 lncRNAs
plot_top50_heatmap <- function(vst_file, deg_file, cancer_type) {
  
  #  Load log-transformed expression data (VST normalized)
  vst_data <- read.csv(vst_file, row.names = 1)
  
  #  Load DE results (shrunken with lfcShrink)
  # Assumes deg_file contains DEGs already filtered for padj < 0.01 and |log2FC| > 1
  degs <- read.csv(deg_file, row.names = 1)
  
  #  Select top 50 genes
  top_genes <- rownames(degs)[1:min(50, nrow(degs))]
  
  #  Subset the expression matrix
  heatmap_data <- vst_data[top_genes, , drop = FALSE]
  
  # Create sample type annotation
  # Labels each sample as either Tumor or Normal based on column names.
  sample_type <- ifelse(grepl("_tumor$", colnames(heatmap_data)), "Tumor", "Normal")
  annotation_col <- data.frame(Sample_Type = sample_type)
  rownames(annotation_col) <- colnames(heatmap_data)
  
  # Reorder columns to group Normal then Tumor
  tumor_samples <- colnames(heatmap_data)[sample_type == "Tumor"]
  normal_samples <- colnames(heatmap_data)[sample_type == "Normal"]
  ordered_cols <- c(normal_samples, tumor_samples)
  
  heatmap_data <- heatmap_data[, ordered_cols]
  annotation_col <- annotation_col[ordered_cols, , drop = FALSE]
  
  # Define annotation colors for the legend
  ann_colors <- list(Sample_Type = c("Tumor" = "#A80000", "Normal" = "#F4A8A8"))
  
  # Save as PNG (high resolution for poster)
  png_filename <- paste0(cancer_type, "_Top50_lncRNA_Heatmap.png")
  png(png_filename, width = 1400, height = 1000, res = 150)
  
  pheatmap(
    mat = heatmap_data,
    scale = "row",                            # Z-score per gene
    show_rownames = TRUE,                     # Show gene names
    show_colnames = FALSE,                    # Hide sample names
    cluster_cols = FALSE,                     # Keep Tumor/Normal order
    cluster_rows = TRUE,                      # Cluster genes
    annotation_col = annotation_col,          # Add Tumor/Normal bar
    annotation_colors = ann_colors,           # Legend colors
    main = paste("Top 50 DE lncRNAs (", cancer_type, ")", sep = ""),
    fontsize = 10,
    fontsize_row = 7,
    color = colorRampPalette(c("#708090", "white", "#A80000"))(100)
  )
  
  dev.off()
}

# LUAD
plot_top50_heatmap(
  vst_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  deg_file = "LUAD_DESeq2_all_results_shrunk.csv",
  cancer_type = "LUAD"
)

# LUSC
plot_top50_heatmap(
  vst_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  deg_file = "LUSC_DESeq2_all_results_shrunk.csv",
  cancer_type = "LUSC"
)


#----------------Venn for LUAD and LUSC lncRNAs---------------------------------
#Load Libraries 
library(dplyr)
library(VennDiagram)


# Load DESeq2 shrunken results (already filtered for p-value & log2FC)
res_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
res_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)


# Extract
sig_LUAD <- rownames(res_LUAD)
sig_LUSC <- rownames(res_LUSC)

#Create the venn
venn <- venn.diagram(
  x = list(
    LUAD = sig_LUAD,
    LUSC = sig_LUSC
  ),
  category.names = c("LUAD DEGs", "LUSC DEGs"),
  filename = "LUAD_LUSC_DEG_Venn.png",
  output = TRUE,
  imagetype = "png",
  
  # Size and resolution
  height = 5500,
  width = 5500,
  resolution = 500,
  
  # Circle appearance
  col = "black",
  fill = c("#FF9999", "#A80000"),
  alpha = 0.6,
  
  # Font sizes and bold styling (no fontfamily)
  cex = 4,                  # Count size
  cat.cex = 4,              # Category label size
  main.cex = 3.9,             # Title size
  fontface = "bold",        # Counts bold
  cat.fontface = "bold",    # Category labels bold
  main.fontface = "bold",   # Title bold
  fontfamily = "sans",         
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  
  # Label positioning
  cat.pos = c(-155, 155),
  cat.dist = c(0.05, 0.05),
  
  # Title
  main = "SHARED AND UNIQUE DE lncRNAs"
)



#-------------------- Select Top 15 DEGs + Map Gene Names + Save CSV ----------------------
#Load libraries
library(biomaRt)

# Load DEGs
deg_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
deg_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)

# Filter significant DEGs first
sig_deg_LUAD <- deg_LUAD %>%
  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  dplyr::arrange(padj, desc(abs(log2FoldChange)))

sig_deg_LUSC <- deg_LUSC %>%
  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  dplyr::arrange(padj, desc(abs(log2FoldChange)))

# Select Top 15 DEGs (Ensembl IDs)
top15_LUAD <- rownames(sig_deg_LUAD)[1:15]
top15_LUSC <- rownames(sig_deg_LUSC)[1:15]

# Save Top 15 DEGs to CSV
write.csv(data.frame(gene_id = top15_LUAD), "LUAD_top15_DEGs.csv", row.names = FALSE)
write.csv(data.frame(gene_id = top15_LUSC), "LUSC_top15_DEGs.csv", row.names = FALSE)

# Print Top 15 for verification
cat("Top 15 LUAD genes:\n")
print(top15_LUAD)

cat("\nTop 15 LUSC genes:\n")
print(top15_LUSC)

# Map Ensembl → gene name
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Prepare all unique Ensembl IDs from LUAD and LUSC
all_ensembl_ids <- unique(c(rownames(deg_LUAD), rownames(deg_LUSC)))
all_ensembl_ids <- gsub("\\..*", "", all_ensembl_ids)  # remove version numbers if present

# Query gene names
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_ensembl_ids,
  mart = mart
)

# Create master gene label vector: Ensembl → HGNC symbol (fallback = Ensembl if missing)
master_gene_label <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)

# Save gene map
write.csv(gene_map, "gene_map_LUAD_LUSC.csv", row.names = FALSE)

# Create mapping vectors for LUAD and LUSC Top 15
LUAD_label_map <- gene_map %>%
  dplyr::filter(ensembl_gene_id %in% gsub("\\..*", "", top15_LUAD))

LUSC_label_map <- gene_map %>%
  dplyr::filter(ensembl_gene_id %in% gsub("\\..*", "", top15_LUSC))


#-----------------------Save ensembl annotation + gene names--------------------
# Load required library
library(readr)

# Load gene map CSV (Ensembl IDs → HGNC symbols)
gene_map <- read_csv("gene_map_LUAD_LUSC.csv", show_col_types = FALSE)

# Print preview to check
head(gene_map)

# Create named vector: Ensembl → HGNC
master_gene_label <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)

# Print a sample of the mapping
cat("Example gene label mapping:\n")
print(head(master_gene_label))

# Save to RData file for future use
save(master_gene_label, file = "master_gene_label.RData")

cat("\n Saved master_gene_label.RData successfully.\n")




#-------------------------------BOXPLOTS OF TOP 15 GENES -----------------------
# Load required libraries
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(biomaRt)

# Load log transformed data
LUAD_vst <- read.csv("LUAD_lncRNA_log_trans_after_norm.csv", row.names = 1)
LUSC_vst <- read.csv("LUSC_lncRNA_log_trans_after_norm.csv", row.names = 1)

# Add gene column
LUAD_vst$Gene <- rownames(LUAD_vst)
LUSC_vst$Gene <- rownames(LUSC_vst)

# Filter to Top 15 genes — use already defined top15_LUAD and top15_LUSC
LUAD_top <- LUAD_vst[LUAD_vst$Gene %in% gsub("\\..*", "", top15_LUAD), ]
LUSC_top <- LUSC_vst[LUSC_vst$Gene %in% gsub("\\..*", "", top15_LUSC), ]

# Melt to long format
LUAD_melt <- melt(LUAD_top, id.vars = "Gene")
LUSC_melt <- melt(LUSC_top, id.vars = "Gene")

# Replace Gene with gene name (fallback = Ensembl ID if missing)
LUAD_melt$Gene <- gsub("\\..*", "", LUAD_melt$Gene)
LUAD_melt$Gene <- ifelse(
  is.na(master_gene_label[LUAD_melt$Gene]) | master_gene_label[LUAD_melt$Gene] == "",
  LUAD_melt$Gene,
  master_gene_label[LUAD_melt$Gene]  
)

LUSC_melt$Gene <- gsub("\\..*", "", LUSC_melt$Gene)
LUSC_melt$Gene <- ifelse(
  is.na(master_gene_label[LUSC_melt$Gene]) | master_gene_label[LUSC_melt$Gene] == "",
  LUSC_melt$Gene,
  master_gene_label[LUSC_melt$Gene]  
)


# Add SampleType column based on name suffix
LUAD_melt$SampleType <- ifelse(grepl("_tumor$", LUAD_melt$variable, ignore.case = TRUE), "Tumor", "Normal")
LUSC_melt$SampleType <- ifelse(grepl("_tumor$", LUSC_melt$variable, ignore.case = TRUE), "Tumor", "Normal")

# Plot LUAD
p_luad <- ggplot(LUAD_melt, aes(x = Gene, y = value, fill = SampleType)) +
  geom_boxplot(
    width = 0.85,           
    outlier.shape = NA,
    alpha = 0.7
  ) +  
  stat_compare_means(
    aes(group = SampleType),
    method = "wilcox.test",
    label = "p.signif",
    size = 12,
    tip.length = 0.01,
    label.y = max(LUAD_melt$value) * 1.2
  ) +
  scale_fill_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
  coord_cartesian(ylim = c(min(LUAD_melt$value), max(LUAD_melt$value) * 1.3)) +
  theme_minimal(base_size = 18) +
  ggtitle("Top 15 Differentially Expressed lncRNAs in LUAD") +
  ylab("VST") +
  xlab("Gene") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  )



ggsave("LUAD_Top15_Boxplot_with_GeneNames.png", plot = p_luad, width = 16, height = 10, dpi = 300)


# Plot LUSC
p_lusc <- ggplot(LUSC_melt, aes(x = Gene, y = value, fill = SampleType)) +
  geom_boxplot(
    width = 0.85,           
    outlier.shape = NA,
    alpha = 0.7
  ) +  
  stat_compare_means(
    aes(group = SampleType),
    method = "wilcox.test",
    label = "p.signif",
    size = 12,
    tip.length = 0.01,
    label.y = max(LUSC_melt$value) * 1.2
  ) +
  scale_fill_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
  coord_cartesian(ylim = c(min(LUSC_melt$value), max(LUSC_melt$value) * 1.3)) +
  theme_minimal(base_size = 18) +
  ggtitle("Top 15 Differentially Expressed lncRNAs in LUSC") +
  ylab("VST") +
  xlab("Gene") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  )

# Save the LUSC plot
ggsave("LUSC_Top15_Boxplot_with_GeneNames.png", plot = p_lusc, width = 16, height = 10, dpi = 300)

# Show both plots
print(p_luad)
print(p_lusc)




#----------------------- Prepare clinical data  --------------------------------
# Load Required Libraries
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(DESeq2)
library(ggplot2)
library(biomaRt)  

# Load normalized counts from DESeq2 output
LUAD_expr <- read.csv("LUAD_lncRNA_normalized_counts_DESeq2.csv", row.names = 1)
LUSC_expr <- read.csv("LUSC_lncRNA_normalized_counts_DESeq2.csv", row.names = 1)

# Load DEGs (for selecting top 15 genes)
deg_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
deg_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)

# Load clinical metadata from TCGA
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD")
clinical_LUSC <- GDCquery_clinic(project = "TCGA-LUSC")

# Convert vital_status into binary outcome
clinical_LUAD$deceased <- ifelse(clinical_LUAD$vital_status == "Alive", FALSE, TRUE)
clinical_LUSC$deceased <- ifelse(clinical_LUSC$vital_status == "Alive", FALSE, TRUE)

# Create overall_survival column
clinical_LUAD$overall_survival <- ifelse(clinical_LUAD$deceased,
                                         clinical_LUAD$days_to_death,
                                         clinical_LUAD$days_to_last_follow_up)

clinical_LUSC$overall_survival <- ifelse(clinical_LUSC$deceased,
                                         clinical_LUSC$days_to_death,
                                         clinical_LUSC$days_to_last_follow_up)

# Fix Expression Sample IDs
fix_sample_ids <- function(names_vec) {
  sapply(names_vec, function(x) {
    parts <- unlist(strsplit(x, split = "\\."))
    parts <- parts[parts != ""] 
    if (length(parts) >= 4) {
      paste(parts[1:3], collapse = "-")
    } else if (length(parts) >= 3) {
      paste("TCGA", parts[2], parts[3], sep = "-")
    } else {
      NA
    }
  })
}

# Apply ID fixing
colnames(LUAD_expr) <- fix_sample_ids(colnames(LUAD_expr))
colnames(LUSC_expr) <- fix_sample_ids(colnames(LUSC_expr))
LUSC_expr <- LUSC_expr[, !is.na(colnames(LUSC_expr))]

# Align expression data and clinical metadata
prepare_survival_data <- function(expr_matrix, clinical_data) {
  shared_samples <- intersect(colnames(expr_matrix), clinical_data$submitter_id)
  
  expr_matrix <- expr_matrix[, shared_samples]
  clinical_data <- clinical_data[clinical_data$submitter_id %in% shared_samples, ]
  clinical_data <- clinical_data[match(colnames(expr_matrix), clinical_data$submitter_id), ]
  
  list(expr = expr_matrix, clinical = clinical_data)
}

LUAD_data <- prepare_survival_data(LUAD_expr, clinical_LUAD)
LUSC_data <- prepare_survival_data(LUSC_expr, clinical_LUSC)

# Save matched data 
LUAD_clinical_clean <- as.data.frame(lapply(LUAD_data$clinical, function(x) {
  if (is.list(x)) sapply(x, toString) else x
}), stringsAsFactors = FALSE)

LUSC_clinical_clean <- as.data.frame(lapply(LUSC_data$clinical, function(x) {
  if (is.list(x)) sapply(x, toString) else x
}), stringsAsFactors = FALSE)

write.csv(LUAD_clinical_clean, "clinical_LUAD_clean.csv", row.names = FALSE)
write.csv(LUSC_clinical_clean, "clinical_LUSC_clean.csv", row.names = FALSE)
write.csv(LUAD_data$expr, "LUAD_expr_matched.csv")
write.csv(LUSC_data$expr, "LUSC_expr_matched.csv")


#--------------------- Match normal samples with clinical data -----------------
# Reload necessary clinical files 
clinical_LUAD <- read.csv("clinical_LUAD_clean.csv")
clinical_LUSC <- read.csv("clinical_LUSC_clean.csv")


# Fix Expression Sample IDs
fix_sample_ids <- function(names_vec) {
  sapply(names_vec, function(x) {
    parts <- unlist(strsplit(x, split = "\\."))
    parts <- parts[parts != ""] 
    if (length(parts) >= 4) {
      paste(parts[1:3], collapse = "-")
    } else if (length(parts) >= 3) {
      paste("TCGA", parts[2], parts[3], sep = "-")
    } else {
      NA
    }
  })
}
# Helper function to check matching samples
check_matching_samples <- function(expr_matrix, clinical_df, label) {
  common_samples <- intersect(colnames(expr_matrix), clinical_df$submitter_id)
  cat("\n==========", label, "==========\n")
  cat("Total expression samples:", ncol(expr_matrix), "\n")
  cat("Matched clinical samples:", length(common_samples), "\n\n")
  return(length(common_samples))
}

# Load Normal expression matrices
LUAD_expr_normal <- read.csv("LUAD_control_lncRNA.csv", row.names = 1)
LUSC_expr_normal <- read.csv("LUSC_control_lncRNA.csv", row.names = 1)

# Fix sample names
colnames(LUAD_expr_normal) <- fix_sample_ids(colnames(LUAD_expr_normal))
colnames(LUSC_expr_normal) <- fix_sample_ids(colnames(LUSC_expr_normal))

# Check sample matching
n_LUAD_normal <- check_matching_samples(LUAD_expr_normal, clinical_LUAD, "LUAD_Normal")
n_LUSC_normal <- check_matching_samples(LUSC_expr_normal, clinical_LUSC, "LUSC_Normal")

# Align expression and clinical data

# LUAD
common_samples_LUAD <- intersect(colnames(LUAD_expr_normal), clinical_LUAD$submitter_id)
LUAD_expr_normal <- LUAD_expr_normal[, common_samples_LUAD]
clinical_LUAD_normal <- clinical_LUAD[match(common_samples_LUAD, clinical_LUAD$submitter_id), ]

# LUSC
common_samples_LUSC <- intersect(colnames(LUSC_expr_normal), clinical_LUSC$submitter_id)
LUSC_expr_normal <- LUSC_expr_normal[, common_samples_LUSC]
clinical_LUSC_normal <- clinical_LUSC[match(common_samples_LUSC, clinical_LUSC$submitter_id), ]

# Save the matched LUAD Normal matrix
write.csv(LUAD_expr_normal, "LUAD_expr_normal_matched.csv")
write.csv(LUSC_expr_normal, "LUSC_expr_normal_matched.csv")

#-----------------Summary for data: number of samples, gender etc---------------
# Load libraries
library(dplyr)
library(readr)
library(knitr)
library(kableExtra)

# Helper to clean expression sample IDs
clean_expr_ids <- function(expr) {
  ids <- gsub("_tumor$|_normal$", "", colnames(expr))
  ids <- gsub("\\.", "-", ids)
  sapply(strsplit(ids, "-"), function(x) paste(x[1:3], collapse = "-"))
}

# Helper to clean clinical data
clean_clinical <- function(df) {
  df$age_years <- df$age_at_diagnosis / 365.25
  df$stage_grouped <- gsub("Stage ([IV]+)[A-B]?", "\\1", df$ajcc_pathologic_stage)
  df$stage_grouped <- factor(df$stage_grouped, levels = c("I", "II", "III", "IV"))
  df$smoking_group <- case_when(
    grepl("(?i)non-smoker|never", df$tobacco_smoking_status) ~ "Never Smoker",
    grepl("(?i)current smoker", df$tobacco_smoking_status) ~ "Current Smoker",
    grepl("(?i)reformed", df$tobacco_smoking_status) ~ "Ex-Smoker",
    TRUE ~ "Unknown"
  )
  return(df)
}

# Load data
expr_LUAD <- read.csv("LUAD_combined.csv", row.names = 1)
expr_LUSC <- read.csv("LUSC_combined.csv", row.names = 1)
clinical_LUAD <- read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)

# Get tumor and normal counts
get_sample_counts <- function(expr) {
  tumor <- sum(grepl("_tumor$", colnames(expr)))
  normal <- sum(grepl("_normal$", colnames(expr)))
  paste0(tumor, " / ", normal)
}
tumor_normal_LUAD <- get_sample_counts(expr_LUAD)
tumor_normal_LUSC <- get_sample_counts(expr_LUSC)

# Match sample IDs
expr_ids_LUAD <- clean_expr_ids(expr_LUAD)
expr_ids_LUSC <- clean_expr_ids(expr_LUSC)
clinical_LUAD$short_id <- substr(clinical_LUAD$submitter_id, 1, 12)
clinical_LUSC$short_id <- substr(clinical_LUSC$submitter_id, 1, 12)

matched_LUAD <- clinical_LUAD[clinical_LUAD$short_id %in% expr_ids_LUAD, ]
matched_LUSC <- clinical_LUSC[clinical_LUSC$short_id %in% expr_ids_LUSC, ]

matched_LUAD <- clean_clinical(matched_LUAD)
matched_LUSC <- clean_clinical(matched_LUSC)

# Build summary table
summary_table <- tibble(
  `Variable` = c(
    "Expression Samples (Tumor / Normal)",
    "Matched Samples (Expression + Clinical)",
    "Age at Diagnosis (mean ± SD)",
    "Age Range (years)",
    "Gender (Male/Female)",
    "Smoking Status (Current/Ex/Never)",
    "Pathologic Stage (I/II/III/IV)",
    "Vital Status (Alive/Dead)"
  ),
  
  `LUAD` = c(
    tumor_normal_LUAD,
    nrow(matched_LUAD),
    sprintf("%.1f ± %.1f", mean(matched_LUAD$age_years, na.rm = TRUE), sd(matched_LUAD$age_years, na.rm = TRUE)),
    sprintf("%.1f – %.1f", min(matched_LUAD$age_years, na.rm = TRUE), max(matched_LUAD$age_years, na.rm = TRUE)),
    paste(table(matched_LUAD$gender), collapse = " / "),
    paste(table(matched_LUAD$smoking_group)[c("Current Smoker", "Ex-Smoker", "Never Smoker")], collapse = " / "),
    paste(table(matched_LUAD$stage_grouped), collapse = " / "),
    paste(table(matched_LUAD$vital_status), collapse = " / ")
  ),
  
  `LUSC` = c(
    tumor_normal_LUSC,
    nrow(matched_LUSC),
    sprintf("%.1f ± %.1f", mean(matched_LUSC$age_years, na.rm = TRUE), sd(matched_LUSC$age_years, na.rm = TRUE)),
    sprintf("%.1f – %.1f", min(matched_LUSC$age_years, na.rm = TRUE), max(matched_LUSC$age_years, na.rm = TRUE)),
    paste(table(matched_LUSC$gender), collapse = " / "),
    paste(table(matched_LUSC$smoking_group)[c("Current Smoker", "Ex-Smoker", "Never Smoker")], collapse = " / "),
    paste(table(matched_LUSC$stage_grouped), collapse = " / "),
    paste(table(matched_LUSC$vital_status), collapse = " / ")
  )
)

# Print table
kbl(summary_table) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))




#----------------Early detection: Boxplots of Top 15 Genes by Stage ------------

library(ggplot2)
library(dplyr)
library(ggpubr)

fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\."))
    parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

plot_stage_vs_normal_top15_facet <- function(expr_tumor_file, expr_normal_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  expr_tumor <- read.csv(expr_tumor_file, row.names = 1)
  expr_normal <- read.csv(expr_normal_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Fix sample IDs
  colnames(expr_tumor) <- fix_sample_ids(colnames(expr_tumor))
  colnames(expr_normal) <- fix_sample_ids(colnames(expr_normal))
  
  # Tumor: group by stage from clinical
  tumor_samples <- colnames(expr_tumor)
  clinical_match <- clinical[match(tumor_samples, clinical$submitter_id), ]
  stage_group <- gsub("^Stage\\s*([IV]+)[A-Z]*$", "Stage \\1", clinical_match$ajcc_pathologic_stage)
  stage_group[is.na(stage_group) | stage_group == ""] <- NA
  # Only keep Stage I and Stage II
  group_tumor <- ifelse(stage_group %in% c("Stage I", "Stage II"), stage_group, NA)
  
  # Normal group
  normal_samples <- colnames(expr_normal)
  group_normal <- rep("Normal", length(normal_samples))
  
  # Restrict to top genes (removing version numbers)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  expr_tumor <- expr_tumor[rownames(expr_tumor) %in% top15_genes_clean, , drop=FALSE]
  expr_normal <- expr_normal[rownames(expr_normal) %in% top15_genes_clean, , drop=FALSE]
  genes_to_plot <- unique(c(rownames(expr_tumor), rownames(expr_normal)))
  
  # Make long-form combined data frame
  all_gene_data <- list()
  for (gene in genes_to_plot) {
    gene_symbol <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    # Tumor
    if (gene %in% rownames(expr_tumor)) {
      expr_vec_tumor <- as.numeric(expr_tumor[gene, ])
      df_tumor <- data.frame(
        expression = expr_vec_tumor,
        group = group_tumor,
        Gene = gene_symbol
      )
      all_gene_data[[paste0(gene, "_tumor")]] <- df_tumor
    }
    # Normal
    if (gene %in% rownames(expr_normal)) {
      expr_vec_normal <- as.numeric(expr_normal[gene, ])
      df_normal <- data.frame(
        expression = expr_vec_normal,
        group = group_normal,
        Gene = gene_symbol
      )
      all_gene_data[[paste0(gene, "_normal")]] <- df_normal
    }
  }
  df <- bind_rows(all_gene_data)
  df <- df[!is.na(df$group), ]
  df$group <- factor(df$group, levels = c("Normal", "Stage I", "Stage II"))
  
  # Sample size per gene/group for n-labels
  label_df <- df %>%
    group_by(Gene, group) %>%
    summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.12, .groups = "drop")
  
  # Define pairwise comparisons for stars
  comparisons_list <- list(
    c("Normal", "Stage I"),
    c("Normal", "Stage II")
  )
  
  # Make plot
  p <- ggplot(df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_text(data = label_df, aes(x = group, y = y, label = paste0("n = ", n)),
              inherit.aes = FALSE, size = 7, color = "black", fontface = "bold") +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.signif",
      size = 9,
      tip.length = 0.03,
      step.increase = 0.2,
      bracket.size = 3,           # Make lines thicker
      bracket.colour = "black"    # Set lines to black for visibility
    ) +scale_fill_manual(values = c(
      "Normal" = "#F4A8A8",
      "Stage I" = "#A0A0A0",
      "Stage II" = "#A80000"
    )) +
    facet_wrap(~Gene, ncol = 4, scales = "free_y") +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("Top 15 lncRNA Expression by Stage (", cancer_type, ")"),
      x = "Group",
      y = "VST"
    ) +
    theme_minimal(base_size = 25) +
    theme(
      plot.title = element_text(size = 48, face = "bold", hjust = 0.5, margin = margin(b = 22)),
      axis.text.x = element_text(angle = 48, hjust = 1, size = 28, face = "bold"),
      axis.title.x = element_text(size = 40, face = "bold"),
      axis.title.y = element_text(size = 44, face = "bold"),
      axis.text.y = element_text(size = 40, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 40, face = "bold"),
      legend.text = element_text(size = 40, face = "bold"),
      strip.text = element_text(size = 40, face = "bold"),
      plot.margin = margin(t = 90, r = 20, b = 20, l = 20)
    )
  
  ggsave(
    paste0("Top15_lncRNAs_by_stage_faceted_", cancer_type, "_EARLY_ONLY_VISIBLE.png"),
    plot = p,
    width = 32,   
    height = 36, 
    dpi = 300
  )
  
  print(p)
  cat("Saved combined faceted plot for", cancer_type, "\n")
}

#for LUAD
plot_stage_vs_normal_top15_facet(
  expr_tumor_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  expr_normal_file = "LUAD_control_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# for LUSC
plot_stage_vs_normal_top15_facet(
  expr_tumor_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  expr_normal_file = "LUSC_control_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)

# ------------------- General Survival Analysis --------------------------------
library(survival)
library(survminer)

# Performs KM survival analysis for each top lncRNA without stratifying by gender/age.

run_general_survival_top15 <- function(expr_matrix, clinical_df, top15_genes, master_gene_label, cancer_type, save_all_plots = TRUE) {
  
  results <- list()
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  for (gene in top15_genes_clean) {
    
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene, ])
    df$deceased <- as.numeric(df$vital_status != "Alive")
    df$overall_survival <- ifelse(df$deceased,
                                  as.numeric(df$days_to_death),
                                  as.numeric(df$days_to_last_follow_up))
    
    df <- df[!is.na(df$expression) & !is.na(df$overall_survival), ]
    
    df$group <- ifelse(df$expression >= median(df$expression, na.rm = TRUE), "High", "Low")
    df <- df[!is.na(df$overall_survival) & !is.na(df$deceased), ]
    
    if (nrow(df) < 10 || length(unique(df$group)) < 2) next
    
    fit <- survfit(Surv(overall_survival, deceased) ~ group, data = df)
    test <- survdiff(Surv(overall_survival, deceased) ~ group, data = df)
    pval <- 1 - pchisq(test$chisq, df = 1)
    
    results[[length(results) + 1]] <- data.frame(
      gene_id = gene,
      gene_symbol = gene_name,
      pval = pval
    )
    
    if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
      plot <- ggsurvplot(
        fit,
        data = df,
        pval = TRUE,
        title = paste(gene_name, "in", cancer_type),
        xlab = "Days",
        ylab = "Survival Probability"
      )
      
      filename <- paste0(cancer_type, "_", gene_name, "_survival.png")
      ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
      
      cat("Saved plot for gene:", gene_name, "\n")
    }
  }
  
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_survival_results_ALL.csv"), row.names = FALSE)
    cat("Saved unstratified survival results for", nrow(results_df), "genes in", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}

run_general_survival_top15(
  expr_matrix = LUAD_expr,
  clinical_df = clinical_LUAD,
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD",
  save_all_plots = TRUE
)

run_general_survival_top15(
  expr_matrix = LUSC_expr,
  clinical_df = clinical_LUSC,
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC",
  save_all_plots = TRUE
)


# --------------------------- Cox Regression for Top 15 Genes ---------------------------
# This function runs univariate Cox proportional hazards regression for each top gene.
# Inputs:
# - expr_matrix: gene expression matrix (genes as rows, samples as columns)
# - clinical_df: clinical dataframe with survival info (vital_status, days_to_death, days_to_last_follow_up)
# - top15_genes: vector of top 15 gene IDs
# - master_gene_label: named vector of gene_symbol names (optional, improves labeling)
# - cancer_type: for naming output file

library(survival)
library(dplyr)

run_grouped_cox_model_top15 <- function(expr_matrix, clinical_df, top15_genes, master_gene_label, cancer_type) {
  
  results <- list()
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  for (gene in top15_genes_clean) {
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene, ])
    df$deceased <- as.numeric(df$vital_status != "Alive")
    df$overall_survival <- ifelse(df$deceased,
                                  as.numeric(df$days_to_death),
                                  as.numeric(df$days_to_last_follow_up))
    
    df <- df[!is.na(df$expression) & !is.na(df$overall_survival), ]
    
    # Create High vs Low expression groups based on median
    df$group <- factor(
      ifelse(df$expression >= median(df$expression, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")  # Set Low as reference
    )
    
    
    if (nrow(df) < 10 || length(unique(df$group)) < 2) next
    
    # Fit Cox model with group as predictor
    cox <- tryCatch({
      summary(coxph(Surv(overall_survival, deceased) ~ group, data = df))
    }, error = function(e) NULL)
    
    if (!is.null(cox)) {
      results[[length(results) + 1]] <- data.frame(
        gene_id = gene,
        gene_symbol = gene_name,
        HR = round(cox$coefficients[,"exp(coef)"], 3),
        lower_CI = round(cox$conf.int[,"lower .95"], 3),
        upper_CI = round(cox$conf.int[,"upper .95"], 3),
        pval = signif(cox$coefficients[,"Pr(>|z|)"], 3)
      )
    }
  }
  
  # Save results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_cox_grouped_model_results.csv"), row.names = FALSE)
    cat("Saved grouped Cox model results for", nrow(results_df), "genes in", cancer_type, "\n")
  } else {
    cat("No valid grouped Cox models were produced.\n")
  }
}

# Read input files first
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1)
clinical_LUAD <- read.csv("clinical_LUAD_clean.csv")

# Now run the function using the data frames (not strings)
run_grouped_cox_model_top15(
  expr_matrix = LUAD_expr,
  clinical_df = clinical_LUAD,
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# Repeat for LUSC
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1)
clinical_LUSC <- read.csv("clinical_LUSC_clean.csv")

run_grouped_cox_model_top15(
  expr_matrix = LUSC_expr,
  clinical_df = clinical_LUSC,
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)


#======Plot cox results
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# Load LUAD and LUSC Cox results
luad <- read_csv("LUAD_cox_grouped_model_results.csv")
lusc <- read_csv("LUSC_cox_grouped_model_results.csv")

# Add cancer type and ensure gene symbols are uppercase
luad <- luad %>% mutate(cancer = "LUAD", gene_symbol = toupper(gene_symbol))
lusc <- lusc %>% mutate(cancer = "LUSC", gene_symbol = toupper(gene_symbol))

# Combine LUAD and LUSC
combined <- bind_rows(luad, lusc)

# Select only your genes of interest
genes_of_interest <- c("BCAN-AS1", "LASTR", "ZFPM2-AS1", "LINC00519")

# Filter to those genes
filtered <- combined %>%
  filter(gene_symbol %in% genes_of_interest) %>%
  mutate(
    label = paste0(gene_symbol, " (", cancer, ")"),
    significance = ifelse(pval < 0.05, "Significant", "Not significant"),
    label = factor(label, levels = rev(label))
  )

# Plot
ggplot(filtered, aes(x = label, y = HR, ymin = lower_CI, ymax = upper_CI)) +
  geom_pointrange(aes(color = significance), size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Significant" = "firebrick", "Not significant" = "grey60")) +
  coord_flip() +
  labs(
    title = "Forest Plot of Selected Prognostic lncRNAs",
    x = "lncRNA (Cancer Type)",
    y = "Hazard Ratio (95% CI)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 20) +  
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )



#------------------------- Correlate Age with Top 15 DEGs ----------------------
# Load required libraries 
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

# Define sample ID fixing function
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Define the correlate_age function (FINAL corrected version)
correlate_age <- function(expr_file, clinical_file, cancer_type, top15_genes, master_gene_label) {
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Fix sample IDs in expression matrix
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match clinical to expression samples
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Convert age from days to years
  clinical$age <- as.numeric(clinical$age_at_diagnosis) / 365
  
  # Run correlation for Top 15 DEGs
  results <- lapply(top15_genes, function(gene) {
    if (gene %in% rownames(expr)) {
      expr_vals <- as.numeric(expr[gene, ])
      age_vals <- clinical$age
      
      if (length(expr_vals) == length(age_vals)) {
        test <- cor.test(expr_vals, age_vals, method = "spearman", exact = FALSE)
        return(data.frame(
          gene_id = gene,
          spearman_rho = test$estimate,
          p_value = test$p.value
        ))
      }
    }
    return(NULL)
  })
  
  # Combine results
  results_df <- do.call(rbind, results)
  
  # If no results, stop safely
  if (is.null(results_df)) {
    cat("No valid correlations computed for", cancer_type, "\n")
    return(NULL)
  }
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Add consistent gene symbols using master_gene_label
  results_df$gene_clean <- gsub("\\..*", "", results_df$gene_id)
  
  results_df$gene_symbol <- ifelse(
    is.na(master_gene_label[results_df$gene_clean]) | master_gene_label[results_df$gene_clean] == "",
    results_df$gene_clean,
    master_gene_label[results_df$gene_clean]
  )
  
  # Save results
  write.csv(results_df, paste0(cancer_type, "_age_correlation_top15.csv"), row.names = FALSE)
  cat("Done with", cancer_type, "- saved age correlation for Top 15 DEGs\n")
}

# Run correlation for Age vs Top 15 DEGs
correlate_age("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label)
correlate_age("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label)

# Age correlation for Normal Samples
correlate_age("LUAD_control_lncRNA.csv", "clinical_LUAD_clean.csv", "LUAD_Normal", top15_LUAD, master_gene_label)
correlate_age("LUSC_control_lncRNA.csv", "clinical_LUSC_clean.csv", "LUSC_Normal", top15_LUSC, master_gene_label)


#==========Plot exp vs age bins

# Load required libraries
library(ggplot2)
library(tidyverse)

#Function for plotting expression vs age bins
plot_age_bin_violin_per_gene_with_significance_flex_v4 <- function(vst_file, clinical_file, cancer_type, top15_genes, master_gene_label, pairwise = FALSE) {
  
  # Load required library for significance stars
  library(ggpubr)
  
  # Load expression and clinical data
  expr <- read.csv(vst_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Clean gene IDs
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Fix sample IDs
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(as.character(x), split = "\\."))
      parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match clinical to expression samples
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Convert age from days to years
  clinical$age <- as.numeric(clinical$age_at_diagnosis) / 365
  
  # For each gene, create 1 violin plot and save it
  for (gene in top15_genes_clean) {
    
    if (gene %in% rownames(expr)) {
      
      expr_vals <- as.numeric(expr[gene, ])
      names(expr_vals) <- NULL  # avoid warning
      
      df <- data.frame(
        expression = expr_vals,
        age = clinical$age  # Use raw numeric age
      )
      
      # Recompute age_bin fresh (no dependency on clinical$age_bin)
      df$age_bin <- cut(
        df$age,
        breaks = c(50, 70, 90),
        right = FALSE,
        labels = c("50-69", "70-89")
      )
      
      # Filter out NA and drop unused factor levels
      df <- df %>% dplyr::filter(!is.na(age_bin))
      df$age_bin <- droplevels(df$age_bin)
      
      # Debug print — see which bins will appear in the plot:
      print(paste0("Bins for ", cancer_type, " - ", gene, ": ", paste(levels(df$age_bin), collapse = ", ")))
      
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Skip empty plots
      if (nrow(df) == 0) {
        message("No valid data for gene ", gene_symbol, " — skipping.")
        next
      }
      
      # Dynamically get available bins
      bins_present <- sort(unique(df$age_bin))
      
      # Create dynamic pairwise comparisons
      if (length(bins_present) >= 2) {
        my_comparisons <- combn(as.character(bins_present), 2, simplify = FALSE)
      } else {
        my_comparisons <- NULL
      }
      
      # Basic plot
      p <- ggplot(df, aes(x = age_bin, y = expression)) +
        geom_violin(trim = FALSE, fill = "#A80000", color = "black", alpha = 0.7) +
        geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5)
      
      # Add significance layer
      if (pairwise && !is.null(my_comparisons)) {
        p <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
          labs(subtitle = "Pairwise Wilcoxon tests")
      } else {
        p <- p + stat_compare_means(method = "kruskal.test", label = "p.signif") +
          labs(subtitle = "Global Kruskal-Wallis test")
      }
      
      # Add titles and theme
      p <- p +
        theme_minimal(base_size = 16) +
        labs(
          title = paste0(cancer_type, " - ", gene_symbol),
          x = "Age Bin (years)",
          y = "Variance Stabilized Expression (VST)"
        ) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14)
        )
      
      # Save plot to PNG
      filename <- if (pairwise && !is.null(my_comparisons)) {
        paste0(cancer_type, "_", gene_symbol, "_AgeBin_PAIRWISE.png")
      } else {
        paste0(cancer_type, "_", gene_symbol, "_AgeBin_GLOBAL.png")
      }
      
      ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
      cat("Saved plot:", filename, "\n")
      
    } else {
      message("Gene not found in expression matrix: ", gene)
    }
  }
  
  cat("Finished saving violin plots for", cancer_type, "\n")
}

#Call the function
plot_age_bin_violin_per_gene_with_significance_flex_v4("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, pairwise = TRUE)
plot_age_bin_violin_per_gene_with_significance_flex_v4("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, pairwise = TRUE)


# Violin plots for Normal
plot_age_bin_violin_per_gene_with_significance_flex_v4("LUAD_control_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD_Normal", top15_LUAD, master_gene_label, pairwise = TRUE)
plot_age_bin_violin_per_gene_with_significance_flex_v4("LUSC_control_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC_Normal", top15_LUSC, master_gene_label, pairwise = TRUE)


#------------------Survival Analysis for DE lncRNAs correlated with age---------
# Load required libraries
library(survival)    # For survival models 
library(survminer)   # For visualization of survival curves with ggplot2


# This function performs survival analysis for each age group (<65 and ≥65 years)
# for each lncRNA significantly correlated with age.
# It saves Kaplan-Meier plots and returns only those with log-rank p < 0.05.
run_age_stratified_survival_all_genes_v2 <- function(expr_matrix, clinical_df, gene_list, master_gene_label, cancer_type, save_all_plots = TRUE) {
  
  # Store all results
  results <- list()
  
  # Loop through all genes
  for (gene in gene_list) {
    gene_clean <- gsub("\\..*", "", gene)
    gene_name <- ifelse(
      is.na(master_gene_label[gene_clean]) | master_gene_label[gene_clean] == "",
      gene_clean,
      master_gene_label[gene_clean]
    )
    
    # Skip if gene not in expression matrix
    if (!(gene_clean %in% rownames(expr_matrix))) next
    
    # Merge expression and clinical
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene_clean, ])
    df$age <- as.numeric(df$age_at_diagnosis) / 365
    df$deceased <- as.numeric(df$vital_status != "Alive")
    df$overall_survival <- ifelse(df$deceased,
                                  as.numeric(df$days_to_death),
                                  as.numeric(df$days_to_last_follow_up))
    
    # Clean missing data
    df <- df[!is.na(df$expression) & !is.na(df$age) & !is.na(df$overall_survival), ]
    
    # Stratify by age group
    df$age_group <- ifelse(df$age < 65, "<65", "≥65")
    
    # Run survival per age group
    for (age_grp in unique(df$age_group)) {
      sub_df <- df[df$age_group == age_grp, ]
      
      # High vs Low expression (median split)
      sub_df$group <- ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low")
      sub_df <- sub_df[!is.na(sub_df$overall_survival) & !is.na(sub_df$deceased), ]
      
      # Skip if too few samples
      if (nrow(sub_df) < 10 || length(unique(sub_df$group)) < 2) next
      
      # Run KM model
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      
      # Log-rank test
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Store result (always)
      results[[length(results)+1]] <- data.frame(
        gene_id = gene_clean,
        gene_symbol = gene_name,
        age_group = age_grp,
        pval = pval
      )
      
      # Save plot — controlled by save_all_plots
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = paste(gene_name, "in", cancer_type, "- Age", age_grp),
          xlab = "Days",
          ylab = "Survival Probability"
        )
        
        ggsave(
          filename = paste0(cancer_type, "_", gene_name, "_Age", gsub("[<>]", "", age_grp), "_survival.png"),
          plot = plot$plot,
          width = 6, height = 5, dpi = 300
        )
      }
    }
  }
  
  # Combine and save all results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_age_stratified_survival_ALL.csv"), row.names = FALSE)
    cat("Saved survival results for", nrow(results_df), "gene/age group combinations for", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}


#Call the function
run_age_stratified_survival_all_genes_v2(LUAD_data$expr, LUAD_data$clinical, top15_LUAD, master_gene_label, "LUAD", save_all_plots = TRUE)
run_age_stratified_survival_all_genes_v2(LUSC_data$expr, LUSC_data$clinical, top15_LUSC, master_gene_label, "LUSC", save_all_plots = TRUE)


# Age-stratified survival for Normal
run_age_stratified_survival_all_genes_v2(LUAD_expr_normal, clinical_LUAD_normal, top15_LUAD, master_gene_label, "LUAD_Normal", save_all_plots = TRUE)
run_age_stratified_survival_all_genes_v2(LUSC_expr_normal, clinical_LUSC_normal, top15_LUSC, master_gene_label, "LUSC_Normal", save_all_plots = TRUE)




#------------------------ Gender-Based Analysis of lncRNAs ------------------------
# Wilcoxon rank-sum test to identify DE lncRNAs significantly different by gender

#Load necessary libraries
library(dplyr)         # Data wrangling
library(readr)         # Reading CSV files
library(AnnotationDbi) # Annotation mapping
library(EnsDb.Hsapiens.v86) # Ensembl annotations
library(rstatix)        # For Wilcoxon test and effect size

# Function to fix sample names
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Combined gender correlation + effect size function
gender_correlation_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Fix sample IDs
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match expression and clinical data
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Clean gene IDs
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Loop over Top 15 genes
  results <- lapply(top15_genes_clean, function(gene) {
    if (gene %in% rownames(expr)) {
      df <- data.frame(
        expression = as.numeric(expr[gene, ]),
        gender = clinical$gender
      ) %>% na.omit()
      
      if (length(unique(df$gender)) == 2) {
        
        # Run Wilcoxon test
        test <- wilcox_test(df, expression ~ gender)
        
        # Run effect size calculation
        effect <- wilcox_effsize(df, expression ~ gender)
        
        # Determine effect direction
        median_male <- median(df$expression[df$gender == "male"], na.rm = TRUE)
        median_female <- median(df$expression[df$gender == "female"], na.rm = TRUE)
        effect_direction <- ifelse(median_male > median_female, "male>female", "female>male")
        
        # Return all results for this gene
        return(tibble(
          gene_id = gene,
          gene_symbol = ifelse(
            is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
            gene,
            master_gene_label[gene]
          ),
          p_value = test$p,
          statistic = test$statistic,
          effsize = effect$effsize,
          magnitude = effect$magnitude,
          effect_direction = effect_direction
        ))
      }
    }
    return(NULL)
  })
  
  # Combine all results
  results_df <- bind_rows(results)
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Save all results
  write.csv(results_df, paste0(cancer_type, "_gender_correlation_Top15.csv"), row.names = FALSE)
  
  cat("Done with", cancer_type, "-", sum(results_df$padj < 0.05), "significant results out of", nrow(results_df), "Top 15 genes tested.\n")
}

# For LUAD
gender_correlation_top15(
  expr_file = "LUAD_expr_matched.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# For LUSC
gender_correlation_top15(
  expr_file = "LUSC_expr_matched.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)

#For normal samples
#LUAD
gender_correlation_top15(
  expr_file = "LUAD_expr_normal_matched.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD_Normal"
)

#LUSC
gender_correlation_top15(
  expr_file = "LUSC_expr_normal_matched.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC_Normal"
)



#====== Plot for Gender-Correlated lncRNAs ======
# Load required libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

#Plotting function
plot_top15_gender_boxplots_for_given_genes <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load log-transformed normalized expression data
  expr <- read.csv(expr_file, row.names = 1)
  
  # Load clinical metadata
  clinical <- read.csv(clinical_file)
  
  # Fix expression sample names to TCGA format
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match samples
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Clean gene IDs in expr
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Use your Top 15 gene list — clean IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Now loop over Top 15 genes
  for (gene in top15_genes_clean) {
    
    if (!(gene %in% rownames(expr))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    expr_vec <- as.numeric(expr[gene, ])
    df <- data.frame(
      expression = expr_vec,
      gender = clinical$gender
    )
    
    # Calculate sample sizes 
    df %>%
      group_by(gender) %>%
      summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.10) -> label_df
    
    label_y_position <- max(df$expression, na.rm = TRUE) * 1.25
    
    # Define gender comparison
    comparisons_list <- list(c("male", "female"))
    
    # Plot
    p <- ggplot(df, aes(x = gender, y = expression, fill = gender)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +  
      stat_compare_means(
        comparisons = comparisons_list,
        method = "wilcox.test",
        label = "p.signif",
        bracket.size = 0.5,
        step.increase = 0.1,  
        size = 5             
      ) +  
      geom_text(
        data = label_df,
        aes(x = gender, y = y, label = paste0("n = ", n)),
        inherit.aes = FALSE,
        size = 5.5,
        color = "black"
      ) +
      scale_fill_manual(values = c("male" = "grey", "female" = "#F8766D")) +
      labs(
        title = paste(gene_name, "expression in", cancer_type),
        x = "Gender",
        y = "Variance Stabilized Expression (VST)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold")
      )
    
    
    ggsave(
      filename = paste0(cancer_type, "_", gene_name, "_gender_boxplot.png"),
      plot = p, width = 8, height = 6, dpi = 300
    )
    
    cat("Saved plot for gene:", gene_name, "\n")
  }
  
  cat("Finished plotting Top 15 gender boxplots for", cancer_type, "\n")
}


# Run the Function for LUAD and LUSC 


# LUAD - Top 15 DEGs
plot_top15_gender_boxplots_for_given_genes(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# LUSC - Top 15 DEGs
plot_top15_gender_boxplots_for_given_genes(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)

# Normal
plot_top15_gender_boxplots_for_given_genes("LUAD_control_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD_Normal")
plot_top15_gender_boxplots_for_given_genes("LUSC_control_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC_Normal")




#--------------------------- Gender-Stratified Survival Analysis ---------------
# This function performs Kaplan-Meier survival analysis separately in male and
# female patients for each lncRNA significantly correlated with gender.

# Required libraries
library(survival)
library(survminer)
library(dplyr)
library(readr)

#Gender-Stratified Survival Function
run_gender_stratified_survival_top15 <- function(expr_matrix, clinical_df, top15_genes, master_gene_label, cancer_type, save_all_plots = TRUE) {
  
  results <- list()  # Store all results
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  for (gene in top15_genes_clean) {
    
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene, ])
    df$deceased <- as.numeric(df$vital_status != "Alive")
    df$overall_survival <- ifelse(df$deceased,
                                  as.numeric(df$days_to_death),
                                  as.numeric(df$days_to_last_follow_up))
    
    df <- df[!is.na(df$expression) & !is.na(df$gender) & !is.na(df$overall_survival), ]
    
    for (sex in unique(df$gender)) {
      
      sub_df <- df[df$gender == sex, ]
      
      sub_df$group <- ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low")
      sub_df <- sub_df[!is.na(sub_df$overall_survival) & !is.na(sub_df$deceased), ]
      
      if (nrow(sub_df) < 10 || length(unique(sub_df$group)) < 2) next
      
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Store result regardless of pval
      results[[length(results) + 1]] <- data.frame(
        gene_id = gene,
        gene_symbol = gene_name,
        gender = sex,
        pval = pval
      )
      
      # Save plot always if save_all_plots == TRUE
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = paste(gene_name, "in", cancer_type, "- Gender", sex),
          xlab = "Days",
          ylab = "Survival Probability"
        )
        
        filename <- paste0(cancer_type, "_", gene_name, "_Gender", sex, "_survival.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        
        cat("Saved plot for gene:", gene_name, "- Gender:", sex, "\n")
      }
    }
  }
  
  # Save results table
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_gender_stratified_survival_ALL.csv"), row.names = FALSE)
    cat("Saved gender-stratified survival results for", nrow(results_df), "gene/gender combinations for", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}

# LUAD
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1)
clinical_LUAD <- read.csv("clinical_LUAD_clean.csv")

run_gender_stratified_survival_top15(
  expr_matrix = LUAD_expr,
  clinical_df = clinical_LUAD,
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD",
  save_all_plots = TRUE
)


# LUSC
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1)
clinical_LUSC <- read.csv("clinical_LUSC_clean.csv")

run_gender_stratified_survival_top15(
  expr_matrix = LUSC_expr,
  clinical_df = clinical_LUSC,
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC",
  save_all_plots = TRUE
)




#------------------------ Smoking-Based Analysis of lncRNAs (3 groups) ------------------------
# Kruskal-Wallis test to identify DE lncRNAs significantly different by smoking status

# Load necessary libraries
library(dplyr)         # Data wrangling
library(readr)         # Reading CSV files
library(AnnotationDbi) # Annotation mapping
library(EnsDb.Hsapiens.v86) # Ensembl annotations

# Function to fix sample names
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Helper function 
add_smoking_group <- function(clinical_df) {
  clinical_df$smoking_group <- dplyr::case_when(
    clinical_df$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
    grepl("Reformed", clinical_df$tobacco_smoking_status) ~ "Ex-smoker",
    clinical_df$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
    TRUE ~ NA_character_
  )
  return(clinical_df)
}

# Function to correlate lncRNA expression with smoking status using Kruskal-Wallis test
smoking_correlation_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Fix expression sample IDs
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Reorder clinical to match expression samples
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Recode smoking status into 3 groups
  clinical <- add_smoking_group(clinical)
  
  # Print counts per group
  print(table(clinical$smoking_group))
  
  # Clean Top 15 gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Apply Kruskal-Wallis test for each gene
  results <- lapply(top15_genes_clean, function(gene) {
    if (gene %in% rownames(expr)) {
      expr_vals <- as.numeric(expr[gene, ])
      smoking_vals <- clinical$smoking_group
      
      if (length(unique(na.omit(smoking_vals))) >= 2) {
        test <- kruskal.test(expr_vals ~ smoking_vals)
        
        group_medians <- tapply(expr_vals, smoking_vals, median, na.rm = TRUE)
        abs_median_diff <- abs(max(group_medians, na.rm = TRUE) - min(group_medians, na.rm = TRUE))
        
        highest_group <- names(group_medians)[which.max(group_medians)]
        lowest_group <- names(group_medians)[which.min(group_medians)]
        direction <- paste0(highest_group, " > ", lowest_group)
        
        return(data.frame(
          gene_id = gene,
          p_value = test$p.value,
          abs_median_diff = abs_median_diff,
          direction = direction
        ))
      }
    }
    return(NULL)
  })
  
  # Combine results
  results_df <- do.call(rbind, results)
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Map gene symbols
  results_df$gene_symbol <- ifelse(
    is.na(master_gene_label[results_df$gene_id]) | master_gene_label[results_df$gene_id] == "",
    results_df$gene_id,
    master_gene_label[results_df$gene_id]
  )
  
  # Save results
  write.csv(results_df, paste0(cancer_type, "_smoking3_correlation_TOP15.csv"), row.names = FALSE)
  cat("Done with", cancer_type, "- saved smoking correlation for Top 15 DEGs\n")
}



# LUAD
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1)
clinical_LUAD <- read.csv("clinical_LUAD_clean.csv")

smoking_correlation_top15(
  expr_file = "LUAD_expr_matched.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)



# LUSC
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1)
clinical_LUSC <- read.csv("clinical_LUSC_clean.csv")

smoking_correlation_top15(
  expr_file = "LUSC_expr_matched.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)


# LUAD Normal
smoking_correlation_top15(
  expr_file = "LUAD_expr_normal_matched.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD_Normal"
)

# LUSC Normal
smoking_correlation_top15(
  expr_file = "LUSC_expr_normal_matched.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC_Normal"
)

#================= Plot Top 15 Smoking-Correlated DEGs =================

# Load required libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)


# Function to fix sample names
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}


# Helper function 
add_smoking_group <- function(clinical_df) {
  clinical_df$smoking_group <- dplyr::case_when(
    clinical_df$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
    grepl("Reformed", clinical_df$tobacco_smoking_status) ~ "Ex-smoker",
    clinical_df$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
    TRUE ~ NA_character_
  )
  return(clinical_df)
}

# Function to plot Top 15 DEGs vs smoking group
plot_top15_smoking_genes <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression data
  expr <- read.csv(expr_file, row.names = 1)
  
  # Load clinical data
  clinical <- read.csv(clinical_file)
  
  # Fix expression sample names
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Recode smoking group
  clinical <- add_smoking_group(clinical)
  
  # Match samples
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Plot for each Top 15 gene
  for (gene in top15_genes_clean) {
    if (gene %in% rownames(expr)) {
      
      expr_vec <- as.numeric(expr[gene, ])
      df <- data.frame(expression = expr_vec, smoking_group = clinical$smoking_group)
      
      # Drop NAs
      df <- df[!is.na(df$smoking_group), ]
      
      # Sample sizes for labels
      df %>%
        group_by(smoking_group) %>%
        summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.10) -> label_df
      
      # Define comparisons
      comparisons_list <- list(
        c("Never smoker", "Ex-smoker"),
        c("Never smoker", "Current smoker"),
        c("Ex-smoker", "Current smoker")
      )
      
      # Adjust significance label positions
      label_y_position <- seq(1.20, 1.10, length.out = length(comparisons_list)) * max(df$expression, na.rm = TRUE)
      
      # Resolve gene name
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Plot
      p <- ggplot(df, aes(x = smoking_group, y = expression, fill = smoking_group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
        stat_compare_means(
          comparisons = comparisons_list,
          method = "wilcox.test",
          label = "p.signif",
          bracket.size = 0.5,
          step.increase = 0.1,
          size = 5
        ) +
        geom_text(data = label_df, aes(x = smoking_group, y = y, label = paste0("n = ", n)),
                  inherit.aes = FALSE, size = 5.5, color = "black") +
        scale_fill_manual(values = c("Never smoker" = "grey", "Ex-smoker" = "#A80000", "Current smoker" = "#F8766D")) +
        coord_cartesian(ylim = c(NA, max(df$expression, na.rm = TRUE) * 1.4)) +  
        labs(
          title = paste(gene_symbol, "expression in", cancer_type),
          x = "Smoking group",
          y = "Variance Stabilized Expression (VST)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 22),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          plot.title = element_text(size = 22, face = "bold")
        )
      
      
      
      # Save plot
      ggsave(
        filename = paste0(cancer_type, "_", gene_symbol, "_Top15_Smoking_Boxplot.png"),
        plot = p, width = 10, height = 8, dpi = 300
      )
      
      cat("Saved plot for gene:", gene_symbol, "\n")
      
    } else {
      message("Gene not found in expression matrix: ", gene)
    }
  }
  
  cat("Finished plotting Top 15 smoking-correlated lncRNAs for", cancer_type, "\n")
}


# LUAD
plot_top15_smoking_genes(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# LUSC 
plot_top15_smoking_genes(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)

  
# Normal
plot_top15_smoking_genes("LUAD_control_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD normal")
plot_top15_smoking_genes("LUSC_control_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC normal")




# ---------------- Smoking-Stratified Survival Analysis - Top 15 ---------------
# Required libraries
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(stringr)

run_smoking_stratified_survival_top15 <- function(expr_matrix, clinical_df, top15_genes, master_gene_label, cancer_type, save_all_plots = TRUE) {
  
  results <- list()  # To store all results
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  for (gene in top15_genes_clean) {
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene, ])
    df$deceased <- as.numeric(df$vital_status != "Alive")
    df$overall_survival <- ifelse(df$deceased,
                                  as.numeric(df$days_to_death),
                                  as.numeric(df$days_to_last_follow_up))
    
    # Recode smoking group
    df$smoking_group <- dplyr::case_when(
      df$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
      grepl("Reformed", df$tobacco_smoking_status) ~ "Ex-smoker",
      df$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
      TRUE ~ NA_character_
    )
    
    # Filter valid data
    df <- df %>%
      dplyr::filter(!is.na(expression), !is.na(smoking_group), !is.na(overall_survival), !is.na(deceased))
    
    for (smoking in unique(df$smoking_group)) {
      sub_df <- df %>% dplyr::filter(smoking_group == smoking)
      
      # High vs Low by median
      sub_df$group <- ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low")
      
      # Filter incomplete
      sub_df <- sub_df %>% dplyr::filter(!is.na(overall_survival), !is.na(deceased), !is.na(group))
    
      
      if (nrow(sub_df) < 10 || length(unique(sub_df$group)) < 2) next
      
      # Survival model
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Save result always
      results[[length(results) + 1]] <- data.frame(
        gene_id = gene,
        gene_symbol = gene_name,
        smoking_group = smoking,
        pval = pval
      )
      
      # Save plot
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        plot_title <- str_wrap(paste0(gene_name, " in ", cancer_type, " - Smoking group: ", smoking), width = 50)
        
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = plot_title,
          xlab = "Days",
          ylab = "Survival Probability"
        )
        
        plot$plot <- plot$plot + theme(
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
          plot.title = element_text(hjust = 0.5, size = 14)
        )
        
        filename <- paste0(cancer_type, "_", gene_name, "_SmokingGroup_", gsub(" ", "_", smoking), "_survival.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        
        cat("Saved plot for gene:", gene_name, "in smoking group:", smoking, "\n")
      }
    }
  }
  
  # Save CSV of all results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_smoking_stratified_survival_ALL.csv"), row.names = FALSE)
    cat("Saved survival results for", nrow(results_df), "gene/smoking group combinations for", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}

# LUAD
run_smoking_stratified_survival_top15(
  LUAD_expr,
  clinical_LUAD,
  top15_LUAD,
  master_gene_label,
  "LUAD",
  save_all_plots = TRUE
)

# LUSC
run_smoking_stratified_survival_top15(
  LUSC_expr,
  clinical_LUSC,
  top15_LUSC,
  master_gene_label,
  "LUSC",
  save_all_plots = TRUE
)




# ------------------------ Stage-Based Analysis of lncRNAs ------------------------
# Kruskal-Wallis test to identify DE lncRNAs significantly different by stage (I, II, III, IV)

# Load necessary libraries
library(dplyr)         # Data wrangling
library(readr)         # Reading CSV files
library(AnnotationDbi) # Annotation mapping
library(EnsDb.Hsapiens.v86) # Ensembl annotations

# Function to fix sample names (for expression matrix)
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Function to perform stage-wise correlation (Kruskal-Wallis test)
stage_correlation_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression data and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Fix sample IDs in expression matrix
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Align expression and clinical samples
  shared_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, shared_samples]
  clinical <- clinical[match(shared_samples, clinical$submitter_id), ]
  
  # Clean ajcc_pathologic_stage → extract Stage I, II, III, IV
  clinical$stage_group <- gsub("^Stage\\s+([IV]+)[A-Z]*$", "Stage \\1", clinical$ajcc_pathologic_stage)
  
  # Print counts of each stage (debug step)
  print(table(clinical$stage_group))
  
  # Clean Top15 gene IDs (remove version if present)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Run Kruskal-Wallis for each gene
  results <- lapply(top15_genes_clean, function(gene) {
    if (gene %in% rownames(expr)) {
      expr_vals <- as.numeric(expr[gene, ])
      stage_vals <- clinical$stage_group
      
      # Only test if at least 2 stages are present
      if (length(unique(na.omit(stage_vals))) >= 2) {
        test <- kruskal.test(expr_vals ~ stage_vals)
        
        # Compute abs median difference
        group_medians <- tapply(expr_vals, stage_vals, median, na.rm = TRUE)
        abs_median_diff <- abs(max(group_medians, na.rm = TRUE) - min(group_medians, na.rm = TRUE))
        
        # Map gene symbol (safe fallback)
        gene_symbol <- ifelse(
          is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
          gene,
          master_gene_label[gene]
        )
        
        return(data.frame(
          gene_id = gene,
          gene_symbol = gene_symbol,
          p_value = test$p.value,
          abs_median_diff = abs_median_diff
        ))
      }
    }
    return(NULL)
  })
  
  # Combine results
  results_df <- do.call(rbind, results)
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    cat("No valid results to process (no variation in stages?)\n")
    return(NULL)
  }
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Save results
  write.csv(results_df, paste0(cancer_type, "_stage_correlation_top15_all.csv"), row.names = FALSE)
  write.csv(results_df[results_df$padj < 0.05, ], paste0(cancer_type, "_stage_correlation_top15_significant.csv"), row.names = FALSE)
  
  cat("Done with", cancer_type, "- Top 15 stage correlation: ", sum(results_df$padj < 0.05), " significant results\n")
}

# LUAD
stage_correlation_top15("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD")

# LUSC
stage_correlation_top15("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC")


#=============Plot stage-wise correlation
# Required libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

# Function to fix sample names
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

plot_top_stage_genes_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression data
  expr <- read.csv(expr_file, row.names = 1)
  
  # Load clinical metadata
  clinical <- read.csv(clinical_file)
  
  # Fix expression sample names to TCGA format
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Clean ajcc_pathologic_stage → Stage I/II/III/IV
  clinical$stage_group <- gsub("^Stage\\s+([IV]+)[A-Z]*$", "Stage \\1", clinical$ajcc_pathologic_stage)
  
  # Match samples
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Clean Top15 gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Plot each Top15 gene
  for (gene in top15_genes_clean) {
    
    if (gene %in% rownames(expr)) {
      
      expr_vec <- as.numeric(as.character(expr[gene, ]))
      df <- data.frame(
        expression = expr_vec,
        stage_group = clinical$stage_group
      )
      
      # Drop NAs
      df <- df[!is.na(df$stage_group), ]
      
      # Skip if < 2 stages available
      if (length(unique(df$stage_group)) < 2) next
      
      # Sample sizes for labels
      df %>%
        group_by(stage_group) %>%
        summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.15) -> label_df
      
      # Define all possible pairwise stage comparisons
      comparisons_list <- list(
        c("Stage I", "Stage II"),
        c("Stage I", "Stage III"),
        c("Stage I", "Stage IV"),
        c("Stage II", "Stage III"),
        c("Stage II", "Stage IV"),
        c("Stage III", "Stage IV")
      )
      
      # Adjust label heights to prevent overlap
      label_positions <- seq(1.40, 1.20, length.out = length(comparisons_list)) * max(df$expression, na.rm = TRUE)
      
      # Map gene symbol
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Plot
      p <- ggplot(df, aes(x = stage_group, y = expression, fill = stage_group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        stat_compare_means(
          comparisons = comparisons_list,
          method = "wilcox.test",
          label = "p.signif",
          tip.length = 0.01,
          step.increase = 0.1,
          size = 5
        ) +
        geom_text(data = label_df, aes(x = stage_group, y = y, label = paste0("n = ", n)),
                  inherit.aes = FALSE, size = 5.5, color = "black") +
        scale_fill_manual(values = c(
          "Stage I" = "red",
          "Stage II" = "grey",
          "Stage III" = "#A80000",
          "Stage IV" = "#D62728"
        )) +
        labs(
          title = paste(gene_symbol, "expression by stage in", cancer_type),
          x = "Pathologic Stage",
          y = "Variance Stabilized Expression (VST)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold")
        )
      
      # Save plot
      ggsave(
        filename = paste0(cancer_type, "_", gene_symbol, "_stage_boxplot_top15.png"),
        plot = p, width = 6, height = 5, dpi = 300
      )
      
      cat("Saved plot for gene:", gene_symbol, "\n")
    } else {
      cat("Gene not found in expression matrix:", gene, "\n")
    }
  }
  
  cat("Finished plotting Top 15 stage-correlated lncRNAs for", cancer_type, "\n")
}


# LUAD
plot_top_stage_genes_top15(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# LUSC
plot_top_stage_genes_top15(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)




#=====================Plotting by grouping stages early vs advanced
# Required libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

# Function to fix sample names
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Function to plot top N stage-correlated lncRNAs
plot_top_stage_genes_grouped_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression data
  expr <- read.csv(expr_file, row.names = 1)
  
  # Load clinical metadata
  clinical <- read.csv(clinical_file)
  
  # Fix expression sample names to TCGA format
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Clean ajcc_pathologic_stage → Stage I/II/III/IV
  clinical$stage_group <- gsub("^Stage\\s+([IV]+)[A-Z]*$", "Stage \\1", clinical$ajcc_pathologic_stage)
  
  # Group stages into Early vs Advanced
  clinical$stage_group <- dplyr::case_when(
    clinical$stage_group %in% c("Stage I", "Stage II") ~ "Early Stage (I/II)",
    clinical$stage_group %in% c("Stage III", "Stage IV") ~ "Advanced Stage (III/IV)",
    TRUE ~ NA_character_
  )
  
  # Match samples
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Clean Top15 gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Plot each Top15 gene
  for (gene in top15_genes_clean) {
    
    if (gene %in% rownames(expr)) {
      
      expr_vec <- as.numeric(as.character(expr[gene, ]))
      df <- data.frame(
        expression = expr_vec,
        stage_group = clinical$stage_group
      )
      
      # Drop NAs
      df <- df[!is.na(df$stage_group), ]
      
      # Skip if < 2 groups available
      if (length(unique(df$stage_group)) < 2) next
      
      # Sample sizes for labels
      df %>%
        group_by(stage_group) %>%
        summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.15) -> label_df
      
      # Define pairwise stage comparison (Early vs Advanced)
      comparisons_list <- list(
        c("Early Stage (I/II)", "Advanced Stage (III/IV)")
      )
      
      # Adjust label heights
      label_positions <- seq(1.30, 1.15, length.out = length(comparisons_list)) * max(df$expression, na.rm = TRUE)
      
      # Map gene symbol
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Plot
      p <- ggplot(df, aes(x = stage_group, y = expression, fill = stage_group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
        geom_jitter(width = 0.2, size = 1.5, alpha = 0.4) +
        stat_compare_means(
          comparisons = comparisons_list,
          method = "wilcox.test",
          label = "p.signif",
          label.y = label_positions
        ) +
        geom_text(data = label_df, aes(x = stage_group, y = y, label = paste0("n = ", n)),
                  inherit.aes = FALSE, size = 4.5, color = "black") +
        scale_fill_manual(values = c(
          "Early Stage (I/II)" = "#A80000",
          "Advanced Stage (III/IV)" = "grey"
        )) +
        labs(
          title = paste(gene_symbol, "expression by stage in", cancer_type),
          x = "Pathologic Stage",
          y = "Log-Transformed Normalized Expression"
        ) +
        theme_minimal(base_size = 14)
      
      # Save plot
      ggsave(
        filename = paste0(cancer_type, "_", gene_symbol, "_stage_grouped_boxplot_top15.png"),
        plot = p, width = 6, height = 5, dpi = 300
      )
      
      cat("Saved plot for gene:", gene_symbol, "\n")
    } else {
      cat("Gene not found in expression matrix:", gene, "\n")
    }
  }
  
  cat("Finished plotting Top 15 stage-grouped lncRNAs for", cancer_type, "\n")
}


# LUAD
plot_top_stage_genes_grouped_top15(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# LUSC
plot_top_stage_genes_grouped_top15(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)



#--------------------Stage-Stratified Survival Analysis for Significant lncRNAs
# Load necessary libraries
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(stringr)

run_stage_stratified_survival_top15 <- function(expr_matrix, clinical_df, top15_genes, master_gene_label, cancer_type, save_all_plots = TRUE) {
  
  results <- list()  # Store all results (for CSV)
  
  # Clean top15 gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  for (gene in top15_genes_clean) {
    
    if (!(gene %in% rownames(expr_matrix))) {
      cat("Gene not found in expression matrix:", gene, "\n")
      next
    }
    
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene, ])
    df$deceased <- as.numeric(df$vital_status != "Alive")
    df$overall_survival <- ifelse(df$deceased,
                                  as.numeric(df$days_to_death),
                                  as.numeric(df$days_to_last_follow_up))
    
    # Clean stage info
    df$stage_group <- gsub("^Stage\\s+([IV]+)[A-Z]*$", "Stage \\1", df$ajcc_pathologic_stage)
    
    # Remove NA rows
    df <- df %>%
      dplyr::filter(!is.na(expression), !is.na(stage_group), !is.na(overall_survival), !is.na(deceased))
    
    # Loop through each stage separately
    for (stage in unique(df$stage_group)) {
      sub_df <- df %>% dplyr::filter(stage_group == stage)
      
      # High/Low expression split (median)
      sub_df$group <- ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low")
      
      # Clean again
      sub_df <- sub_df %>%
        dplyr::filter(!is.na(overall_survival), !is.na(deceased), !is.na(group))
      
      # Skip too small groups or one group only
      if (nrow(sub_df) < 10 || length(unique(sub_df$group)) < 2) next
      
      # Survival model
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Map gene symbol
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Save plot (always if save_all_plots = TRUE or if p < 0.05)
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        
        plot_title <- str_wrap(paste0(gene_symbol, " in ", cancer_type, " - Stage: ", stage), width = 50)
        
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = plot_title,
          xlab = "Days",
          ylab = "Survival Probability"
        )
        
        plot$plot <- plot$plot + theme(
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
          plot.title = element_text(hjust = 0.5, size = 14)
        )
        
        # Save PNG
        filename <- paste0(cancer_type, "_", gene_symbol, "_Stage_", gsub(" ", "_", stage), "_survival_top15.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        
        cat("Saved plot:", filename, "\n")
      }
      
      # Save results row
      results[[length(results) + 1]] <- data.frame(
        gene_id = gene,
        gene_symbol = gene_symbol,
        stage_group = stage,
        pval = pval
      )
    }
  }
  
  # Save summary CSV
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_stage_stratified_survival_ALL_top15.csv"), row.names = FALSE)
    cat("Saved survival results for", nrow(results_df), "stage groups in", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}


# LUAD
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1)
clinical_LUAD <- read.csv("clinical_LUAD_clean.csv")

run_stage_stratified_survival_top15(
  LUAD_expr,
  clinical_LUAD,
  top15_LUAD,
  master_gene_label,
  "LUAD",
  save_all_plots = TRUE
)

# LUSC
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1)
clinical_LUSC <- read.csv("clinical_LUSC_clean.csv")

run_stage_stratified_survival_top15(
  LUSC_expr,
  clinical_LUSC,
  top15_LUSC,
  master_gene_label,
  "LUSC",
  save_all_plots = TRUE
)



#---------------------------- Multivariate cox regression ----------------------
# Load required libraries
library(survival)
library(dplyr)
library(readr)
library(stringr)

# Fix sample IDs
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(as.character(x), split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Multivariate Cox regression for top 15 lncRNAs
run_multivariate_cox <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Fix expression sample IDs
  colnames(expr) <- fix_sample_ids(colnames(expr))
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Match samples
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Prepare clinical variables
  clinical$age <- as.numeric(clinical$age_at_diagnosis) / 365
  clinical$gender <- as.factor(clinical$gender)
  clinical$smoking_group <- dplyr::case_when(
    clinical$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
    grepl("Reformed", clinical$tobacco_smoking_status) ~ "Ex-smoker",
    clinical$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
    TRUE ~ NA_character_
  )
  clinical$smoking_group <- as.factor(clinical$smoking_group)
  clinical$stage_group <- gsub("^Stage\\s+([IV]+)[A-Z]*$", "Stage \\1", clinical$ajcc_pathologic_stage)
  clinical$stage_group <- as.factor(clinical$stage_group)
  clinical$deceased <- as.numeric(clinical$vital_status != "Alive")
  clinical$overall_survival <- ifelse(clinical$deceased,
                                      as.numeric(clinical$days_to_death),
                                      as.numeric(clinical$days_to_last_follow_up))
  
  # Initialize results
  results_list <- list()
  
  for (gene in top15_genes_clean) {
    if (!(gene %in% rownames(expr))) {
      cat("Gene not found:", gene, "\n")
      next
    }
    
    gene_symbol <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    clinical$expression <- as.numeric(expr[gene, ])
    
    df <- clinical %>%
      dplyr::filter(!is.na(expression), !is.na(age), !is.na(gender), 
                    !is.na(smoking_group), !is.na(stage_group), 
                    !is.na(overall_survival), !is.na(deceased))
    
    if (nrow(df) < 20) {
      cat("Too few samples for:", gene_symbol, "\n")
      next
    }
    
    # Fit multivariate Cox model
    cox_model <- tryCatch({
      coxph(Surv(overall_survival, deceased) ~ expression + age + gender + smoking_group + stage_group, data = df)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(cox_model)) {
      summary_cox <- summary(cox_model)
      coefs <- as.data.frame(summary_cox$coefficients)
      coefs$gene_id <- gene
      coefs$gene_symbol <- gene_symbol
      coefs$variable <- rownames(coefs)
      rownames(coefs) <- NULL
      results_list[[gene]] <- coefs
    }
  }
  
  # Combine and save results
  if (length(results_list) > 0) {
    all_results <- bind_rows(results_list)
    write.csv(all_results, paste0(cancer_type, "_multivariate_cox_results_top15.csv"), row.names = FALSE)
    cat("Saved multivariate Cox results for", cancer_type, "\n")
  } else {
    cat("No valid multivariate Cox models for", cancer_type, "\n")
  }
}


# Run for LUAD
run_multivariate_cox("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD")

# Run for LUSC
run_multivariate_cox("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC")


#================= Plot multivariate cox (Ordered by Cancer Type) =================
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)

# Load multivariate Cox results
df_luad <- read_csv("LUAD_multivariate_cox_results_top15.csv")
df_lusc <- read_csv("LUSC_multivariate_cox_results_top15.csv")

# Add cancer type
df_luad$cancer <- "LUAD"
df_lusc$cancer <- "LUSC"

# Combine and keep only 'expression' rows
combined_df <- bind_rows(df_luad, df_lusc) %>%
  filter(variable == "expression")

# Prepare plotting data
plot_df <- combined_df %>%
  mutate(
    HR = exp(coef),
    lower_CI = exp(coef - 1.96 * `se(coef)`),
    upper_CI = exp(coef + 1.96 * `se(coef)`),
    gene_display = paste0(gene_symbol, " (", cancer, ")"),
    pval_label = ifelse(`Pr(>|z|)` < 0.001, "<0.001", sprintf("p = %.3f", `Pr(>|z|)`)),
    cancer = factor(cancer, levels = c("LUAD", "LUSC")) # Force LUAD first
  ) %>%
  arrange(cancer, desc(HR)) %>%  # Sort by cancer then HR
  mutate(gene_display = fct_inorder(gene_display))

# Plot
ggplot(plot_df, aes(x = HR, y = gene_display)) +
  geom_point(size = 5, color = "#A80000") +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10(
    limits = c(
      min(plot_df$lower_CI, na.rm = TRUE),
      max(plot_df$upper_CI, na.rm = TRUE) + 0.001
    ),
    expand = expansion(mult = c(0.01, 0.2))
  ) +
  geom_text(
    aes(label = pval_label, x = upper_CI + 0.0002),
    size = 5,
    hjust = 0
  ) +
  labs(
    title = "Multivariate Cox Model",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  )

# Save as portrait
ggsave("multivariate_forest_plot.png", width = 8, height = 12, dpi = 300)





#------------------Combinations of lncRNAs for predictive potential-------------

# Load required packages
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(combinat)

# Helper: fix sample IDs
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(as.character(x), split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Core function to run combinations
run_combinations <- function(expr_file, clinical_file, gene_list, master_gene_label, cancer_type) {
  
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  colnames(expr) <- fix_sample_ids(colnames(expr))
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  gene_list <- gene_list[gene_list %in% rownames(expr)]  # Ensure valid genes
  
  # Match samples
  shared <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, shared]
  clinical <- clinical[match(shared, clinical$submitter_id), ]
  
  # Preprocess clinical
  clinical$age <- as.numeric(clinical$age_at_diagnosis) / 365
  clinical$gender <- as.factor(clinical$gender)
  clinical$smoking_group <- dplyr::case_when(
    clinical$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
    grepl("Reformed", clinical$tobacco_smoking_status) ~ "Ex-smoker",
    clinical$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
    TRUE ~ NA_character_
  )
  clinical$smoking_group <- as.factor(clinical$smoking_group)
  clinical$stage_group <- gsub("^Stage\\s+([IV]+)[A-Z]*$", "Stage \\1", clinical$ajcc_pathologic_stage)
  clinical$stage_group <- as.factor(clinical$stage_group)
  clinical$deceased <- as.numeric(clinical$vital_status != "Alive")
  clinical$overall_survival <- ifelse(clinical$deceased,
                                      as.numeric(clinical$days_to_death),
                                      as.numeric(clinical$days_to_last_follow_up))
  
  results <- list()
  
  for (k in 2:3) {
    combos <- combn(gene_list, k, simplify = FALSE)
    
    for (genes in combos) {
      df <- clinical
      for (gene in genes) {
        df[[gene]] <- as.numeric(expr[gene, ])
      }
      
      # Filter valid samples
      df <- df[complete.cases(df[, genes]), ]
      df <- df %>%
        filter(!is.na(overall_survival), !is.na(deceased), !is.na(age),
               !is.na(gender), !is.na(smoking_group), !is.na(stage_group))
      
      if (nrow(df) < 20) next
      
      gene_formula <- paste(genes, collapse = " + ")
      model <- tryCatch({
        coxph(as.formula(paste("Surv(overall_survival, deceased) ~", gene_formula)), data = df)
      }, error = function(e) return(NULL))
      
      if (!is.null(model)) {
        p_val <- summary(model)$coefficients[, "Pr(>|z|)"]
        min_p <- min(p_val, na.rm = TRUE)
        if (min_p < 0.05) {
          symbol_names <- sapply(genes, function(x) {
            sym <- master_gene_label[[x]]
            ifelse(is.na(sym) || sym == "", x, sym)
          })
          
          results[[length(results) + 1]] <- list(
            ensembl_genes = genes,
            symbols = symbol_names,
            min_pvalue = min_p
          )
        }
      }
    }
  }
  
  if (length(results) > 0) {
    df_results <- tibble::tibble(
      ensembl_combination = sapply(results, function(x) paste(x$ensembl_genes, collapse = " + ")),
      symbol_combination = sapply(results, function(x) paste(x$symbols, collapse = " + ")),
      min_pvalue = sapply(results, function(x) x$min_pvalue)
    )
    
    out_file <- paste0("significant_", cancer_type, "_lncRNA_combinations.csv")
    write.csv(df_results, out_file, row.names = FALSE)
    cat("Significant combinations saved for", cancer_type, "\n")
  } else {
    cat("No significant combinations found for", cancer_type, "\n")
  }
}

# LUAD
run_combinations("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD")

# LUSC
run_combinations("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC")




#======================Plot results
# Required libraries
library(survival)
library(survminer)
library(dplyr)
library(readr)

# Fix sample IDs
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    parts <- unlist(strsplit(as.character(x), split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Function to plot top 3 KM curves based on combination results
plot_top_combinations <- function(expr_file, clinical_file, combo_results_file, cancer_type, output_dir = ".") {
  # Load data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  combo <- read.csv(combo_results_file)
  
  # Sample ID cleanup
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Keep top 3 based on min_pvalue
  combo <- combo[order(combo$min_pvalue), ]
  top3 <- head(combo, 20)
  
  for (i in 1:nrow(top3)) {
    genes <- unlist(strsplit(top3$ensembl_combination[i], split = " \\+ "))
    genes <- genes[genes %in% rownames(expr)]
    
    if (length(genes) < 2) next  # skip invalid
    
    shared <- intersect(colnames(expr), clinical$submitter_id)
    expr_sub <- expr[genes, shared]
    clinical_sub <- clinical[match(shared, clinical$submitter_id), ]
    
    # Survival info
    clinical_sub$age <- as.numeric(clinical_sub$age_at_diagnosis) / 365
    clinical_sub$deceased <- as.numeric(clinical_sub$vital_status != "Alive")
    clinical_sub$overall_survival <- ifelse(
      clinical_sub$deceased,
      as.numeric(clinical_sub$days_to_death),
      as.numeric(clinical_sub$days_to_last_follow_up)
    )
    
    # Add expression values
    for (gene in genes) {
      clinical_sub[[gene]] <- as.numeric(expr_sub[gene, ])
    }
    
    df <- clinical_sub %>%
      filter(!is.na(overall_survival), !is.na(deceased)) %>%
      filter(rowSums(is.na(across(all_of(genes)))) == 0)
    
    if (nrow(df) < 20) next
    
    # Risk score
    model_formula <- as.formula(paste("Surv(overall_survival, deceased) ~", paste(genes, collapse = " + ")))
    cox_model <- coxph(model_formula, data = df)
    df$risk_score <- scale(predict(cox_model, type = "lp"))
    df$risk_group <- ifelse(df$risk_score >= median(df$risk_score, na.rm = TRUE), "High", "Low")
    
    # Plot
    surv_fit <- survfit(Surv(overall_survival, deceased) ~ risk_group, data = df)
    symbol_label <- gsub(" \\+ ", "_", top3$symbol_combination[i])
    plot_title <- paste0("KM - ", cancer_type, ": ", top3$symbol_combination[i])
    
    km_plot <- ggsurvplot(
      surv_fit,
      data = df,
      pval = TRUE,
      risk.table = TRUE,
      title = plot_title,
      xlab = "Days",
      ylab = "Survival Probability",
      ggtheme = theme_minimal(base_size = 18),
      font.title = c(22, "bold"),
      font.x = c(20),
      font.y = c(20),
      font.tickslab = c(18),
      font.legend = c(18),
      font.subtitle = c(20, "plain"),
      font.caption = c(18),
      risk.table.fontsize = 5
    )
    
    
    filename <- paste0(output_dir, "/", cancer_type, "_KM_top", i, "_", symbol_label, ".png")
    ggsave(filename, plot = km_plot$plot, width = 10, height = 6, dpi = 300)
    cat("Saved plot:", filename, "\n")
  }
}

#Run
plot_top_combinations("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", "significant_LUAD_lncRNA_combinations.csv", "LUAD")
plot_top_combinations("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", "significant_LUSC_lncRNA_combinations.csv", "LUSC")




