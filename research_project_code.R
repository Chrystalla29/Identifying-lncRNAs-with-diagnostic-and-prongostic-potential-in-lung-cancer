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
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
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
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

#Save the plot
ggsave("MA_plot_LUSC.png", plot = p2, width = 8, height = 6, dpi = 300)


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

#Create Venn Diagram 
venn_out <- venn.diagram(
  x = list(
    LUAD = sig_LUAD,
    LUSC = sig_LUSC
  ),
  category.names = c("LUAD DEGs", "LUSC DEGs"),
  filename = "LUAD_LUSC_DEG_Venn.png",  # Save as PNG
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 500,
  col = "black",
  fill = c("#FF9999", "#A80000"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  main = "Shared and Unique DE lncRNAs",
  main.cex = 2
)


#-------------------------------BOXPLOTS OF TOP 5 GENES-------------------------
# Load required libraries
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Load log transformed data and DEGs
LUAD_vst <- read.csv("LUAD_lncRNA_log_trans_after_norm.csv", row.names = 1)
LUSC_vst <- read.csv("LUSC_lncRNA_log_trans_after_norm.csv", row.names = 1)
deg_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
deg_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)

# Sort by lowest padj and highest absolute log2 fold change
deg_LUAD <- deg_LUAD %>% arrange(padj, desc(abs(log2FoldChange)))
deg_LUSC <- deg_LUSC %>% arrange(padj, desc(abs(log2FoldChange)))

#Select Top 5 DEGs
top5_LUAD <- rownames(deg_LUAD)[1:5]
top5_LUSC <- rownames(deg_LUSC)[1:5]

print(top5_LUAD)
print(top5_LUSC)

# Reshape for ggplot
# Add gene column
LUAD_vst$Gene <- rownames(LUAD_vst)
LUSC_vst$Gene <- rownames(LUSC_vst)

# Filter to top genes
LUAD_top <- LUAD_vst[LUAD_vst$Gene %in% top5_LUAD, ]
LUSC_top <- LUSC_vst[LUSC_vst$Gene %in% top5_LUSC, ]

# Melt to long format
LUAD_melt <- melt(LUAD_top, id.vars = "Gene")
LUSC_melt <- melt(LUSC_top, id.vars = "Gene")

# Add SampleType column based on name suffix
LUAD_melt$SampleType <- ifelse(grepl("_tumor$", LUAD_melt$variable, ignore.case = TRUE), "Tumor", "Normal")
LUSC_melt$SampleType <- ifelse(grepl("_tumor$", LUSC_melt$variable, ignore.case = TRUE), "Tumor", "Normal")

#  Plot LUAD 
p_luad <- ggplot(LUAD_melt, aes(x = Gene, y = value, fill = SampleType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(aes(group = SampleType), method = "wilcox.test", label = "p.signif") +  # ⭐ add this line
  scale_fill_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
  theme_minimal(base_size = 13) +
  ggtitle("Top 5 DEGs in LUAD") +
  ylab("Log Transformed Count") +
  xlab("Gene") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

ggsave("LUAD_Top5_Boxplot.png", plot = p_luad, width = 8, height = 5, dpi = 300)

# Plot LUSC 
p_lusc <- ggplot(LUSC_melt, aes(x = Gene, y = value, fill = SampleType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(aes(group = SampleType), method = "wilcox.test", label = "p.signif") +  # ⭐ add here too
  scale_fill_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
  theme_minimal(base_size = 13) +
  ggtitle("Top 5 DEGs in LUSC") +
  ylab("Log Transformed Count") +
  xlab("Gene") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )


ggsave("LUSC_Top5_Boxplot.png", plot = p_lusc, width = 8, height = 5, dpi = 300)

# Show both plots
print(p_luad)
print(p_lusc)




#----------------------Survival analysis for top 5 genes------------------------
# Load Required Libraries for survival analysis and RNA-seq processing
library(TCGAbiolinks)   # Access TCGA data
library(survival)       # Survival analysis
library(survminer)      # Survival plots
library(dplyr)          # Data manipulation
library(readr)          # CSV reading
library(DESeq2)         # Normalization and differential expression
library(ggplot2)        # Plotting

# Load normalized counts from DESeq2 output
LUAD_expr <- read.csv("LUAD_lncRNA_normalized_counts_DESeq2.csv", row.names = 1)
LUSC_expr <- read.csv("LUSC_lncRNA_normalized_counts_DESeq2.csv", row.names = 1)

# Load clinical metadata from TCGA
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD")
clinical_LUSC <- GDCquery_clinic(project = "TCGA-LUSC")

# Convert vital_status into binary outcome: TRUE = dead, FALSE = alive
clinical_LUAD$deceased <- ifelse(clinical_LUAD$vital_status == "Alive", FALSE, TRUE)
clinical_LUSC$deceased <- ifelse(clinical_LUSC$vital_status == "Alive", FALSE, TRUE)

# Create overall_survival column based on death or follow-up time
# Create the time to event variable
clinical_LUAD$overall_survival <- ifelse(clinical_LUAD$deceased,
                                         clinical_LUAD$days_to_death,
                                         clinical_LUAD$days_to_last_follow_up)

clinical_LUSC$overall_survival <- ifelse(clinical_LUSC$deceased,
                                         clinical_LUSC$days_to_death,
                                         clinical_LUSC$days_to_last_follow_up)

# Fix Expression Sample IDs to match clinical submitter_id
fix_sample_ids <- function(names_vec) {
  sapply(names_vec, function(x) {
    parts <- unlist(strsplit(x, split = "\\."))
    parts <- parts[parts != ""]  # Remove empty strings
    if (length(parts) >= 4) {
      paste(parts[1:3], collapse = "-")  # TCGA-XX-XXXX format
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
LUSC_expr <- LUSC_expr[, !is.na(colnames(LUSC_expr))]  # Remove problematic columns

# Align expression data and clinical metadata by common sample IDs
prepare_survival_data <- function(expr_matrix, clinical_data) {
  shared_samples <- intersect(colnames(expr_matrix), clinical_data$submitter_id)
  
  expr_matrix <- expr_matrix[, shared_samples]
  clinical_data <- clinical_data[clinical_data$submitter_id %in% shared_samples, ]
  clinical_data <- clinical_data[match(colnames(expr_matrix), clinical_data$submitter_id), ]
  
  list(expr = expr_matrix, clinical = clinical_data)
}

# Apply matching for LUAD and LUSC datasets
LUAD_data <- prepare_survival_data(LUAD_expr, clinical_LUAD)
LUSC_data <- prepare_survival_data(LUSC_expr, clinical_LUSC)

# Define top 5 most significant DEGs to test survival on
top5_LUAD <- c("ENSG00000227066", "ENSG00000262772", "ENSG00000251003",
               "ENSG00000243479", "ENSG00000227036")

top5_LUSC <- c("ENSG00000249395", "ENSG00000230937", "ENSG00000225548", 
               "ENSG00000224271", "ENSG00000267284")

# Survival Analysis Function
run_survival_analysis <- function(gene, expr_data, clinical_data, cancer_type) {
  
  # Check if gene exists in matrix
  if (!(gene %in% rownames(expr_data))) {
    message("Gene not found in expression matrix: ", gene)
    return(NULL)
  }
  
  # Extract expression values and merge with clinical
  expr_vector <- as.numeric(expr_data[gene, ])
  df <- clinical_data
  
  # Ensure sample sizes match
  if (length(expr_vector) != nrow(df)) {
    message("Sample size mismatch for gene ", gene)
    return(NULL)
  }
  
  df$expression <- expr_vector
  
  # Remove missing data
  df <- df[!is.na(df$expression) & !is.na(df$overall_survival) & !is.na(df$deceased), ]
  
  # Define high and low expression using quantile split (75th percentile)
  high_threshold <- quantile(df$expression, 0.75)
  low_threshold <- quantile(df$expression, 0.25)
  
  # Filter for only extreme high and low expression samples
  df <- df[df$expression >= high_threshold | df$expression <= low_threshold, ]
  
  # Assign groups for survival analysis
  df$group <- ifelse(df$expression >= high_threshold, "High", "Low")
  
  # Skip if only one group remains
  if (length(unique(df$group)) < 2) {
    message("Only one expression group for gene ", gene)
    return(NULL)
  }
  
  # Ensure correct format
  df$group <- as.factor(df$group)
  df$deceased <- as.numeric(df$deceased)
  df$overall_survival <- as.numeric(df$overall_survival)
  
  # Fit KM survival curve and log-rank test
  tryCatch({
    fit <- survfit(Surv(overall_survival, deceased) ~ group, data = df)
    surv_plot <- ggsurvplot(
      fit, data = df,
      pval = TRUE,
      title = paste("Survival for", gene, "in", cancer_type),
      xlab = "Days", ylab = "Survival Probability"
    )
    
    ggsave(
      filename = paste0(cancer_type, "_", gene, "_survival_quantile.png"),
      plot = surv_plot$plot,
      width = 6, height = 5, dpi = 300
    )
    cat("Plot saved for gene:", gene, "\n")
    
  }, error = function(e) {
    message("Error during plotting for gene ", gene, ": ", e$message)
  })
}

# Run Survival for LUAD
for (gene in top5_LUAD) {
  run_survival_analysis(gene, LUAD_data$expr, LUAD_data$clinical, "LUAD")
}

# Run Survival for LUSC
for (gene in top5_LUSC) {
  run_survival_analysis(gene, LUSC_data$expr, LUSC_data$clinical, "LUSC")
}


#--------------------Survival analysis on all DE lncRNAs------------------------
# Load Required Libraries
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(DESeq2)
library(ggplot2)

# Load normalized counts
LUAD_expr <- read.csv("LUAD_lncRNA_normalized_counts_DESeq2.csv", row.names = 1)
LUSC_expr <- read.csv("LUSC_lncRNA_normalized_counts_DESeq2.csv", row.names = 1)

# Fix sample IDs to match TCGA clinical sample format
fix_sample_ids <- function(names_vec) {
  sapply(names_vec, function(x) {
    parts <- unlist(strsplit(x, split = "\\."))
    parts <- parts[parts != ""]
    if (length(parts) >= 4) paste(parts[1:3], collapse = "-")
    else if (length(parts) >= 3) paste("TCGA", parts[2], parts[3], sep = "-")
    else NA
  })
}

# Apply sample ID fixer
colnames(LUAD_expr) <- fix_sample_ids(colnames(LUAD_expr))
colnames(LUSC_expr) <- fix_sample_ids(colnames(LUSC_expr))
LUAD_expr <- LUAD_expr[, !is.na(colnames(LUAD_expr))]
LUSC_expr <- LUSC_expr[, !is.na(colnames(LUSC_expr))]

# Load clinical data
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD")
clinical_LUSC <- GDCquery_clinic(project = "TCGA-LUSC")

# Create binary survival and time columns
clinical_LUAD$deceased <- clinical_LUAD$vital_status != "Alive"
clinical_LUSC$deceased <- clinical_LUSC$vital_status != "Alive"
clinical_LUAD$overall_survival <- ifelse(clinical_LUAD$deceased,
                                         clinical_LUAD$days_to_death,
                                         clinical_LUAD$days_to_last_follow_up)
clinical_LUSC$overall_survival <- ifelse(clinical_LUSC$deceased,
                                         clinical_LUSC$days_to_death,
                                         clinical_LUSC$days_to_last_follow_up)

# Match expression and clinical data
prepare_survival_data <- function(expr_matrix, clinical_data) {
  shared_samples <- intersect(colnames(expr_matrix), clinical_data$submitter_id)
  expr_matrix <- expr_matrix[, shared_samples]
  clinical_data <- clinical_data[clinical_data$submitter_id %in% shared_samples, ]
  clinical_data <- clinical_data[match(colnames(expr_matrix), clinical_data$submitter_id), ]
  list(expr = expr_matrix, clinical = clinical_data)
}

# Prepare LUAD and LUSC datasets
LUAD_data <- prepare_survival_data(LUAD_expr, clinical_LUAD)
LUSC_data <- prepare_survival_data(LUSC_expr, clinical_LUSC)

# Load DE lncRNAs 
res_LUAD_filtered <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
sig_LUAD_genes <- rownames(res_LUAD_filtered)
res_LUSC_filtered <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)
sig_LUSC_genes <- rownames(res_LUSC_filtered)


# Survival analysis function with quantile-based grouping
run_survival_analysis <- function(gene, expr_data, clinical_data, cancer_type) {
  if (!(gene %in% rownames(expr_data))) return(NULL)
  
  expr_vector <- as.numeric(expr_data[gene, ])
  df <- as.data.frame(clinical_data)
  
  if (length(expr_vector) != nrow(df)) return(NULL)
  
  df$expression <- expr_vector
  df <- df[!is.na(df$expression) & !is.na(df$overall_survival) & !is.na(df$deceased), ]
  if (nrow(df) < 5) return(NULL)
  
  # Use quantile split
  q1 <- quantile(df$expression, 0.25)
  q3 <- quantile(df$expression, 0.75)
  df <- df[df$expression <= q1 | df$expression >= q3, ]
  df$group <- ifelse(df$expression >= q3, "High", "Low")
  if (length(unique(df$group)) < 2) return(NULL)
  
  df$group <- as.factor(df$group)
  df$deceased <- as.numeric(df$deceased)
  df$overall_survival <- as.numeric(df$overall_survival)
  
  fit <- tryCatch(survfit(Surv(overall_survival, deceased) ~ group, data = df), error = function(e) return(NULL))
  if (is.null(fit)) return(NULL)
  
  test <- tryCatch(survdiff(Surv(overall_survival, deceased) ~ group, data = df), error = function(e) return(NULL))
  if (is.null(test)) return(NULL)
  
  pval <- 1 - pchisq(test$chisq, df = 1)
  
  # Determine survival direction using named group medians
  median_vals <- summary(fit)$table[, "median"]
  high_group <- which(names(median_vals) == "group=High")
  low_group <- which(names(median_vals) == "group=Low")
  
  if (length(high_group) == 1 && length(low_group) == 1) {
    if (!is.na(median_vals[high_group]) && !is.na(median_vals[low_group])) {
      survival_direction <- ifelse(median_vals[high_group] > median_vals[low_group],
                                   "High expression = better survival",
                                   "High expression = worse survival")
    } else {
      survival_direction <- NA
    }
  } else {
    survival_direction <- NA
  }
  
  # Save and return only if significant
  if (!is.na(pval) && pval < 0.05) {
    surv_plot <- ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      title = paste("Survival for", gene, "in", cancer_type),
      xlab = "Days", ylab = "Survival Probability"
    )
    ggsave(
      filename = paste0(cancer_type, "_", gene, "_survival.png"),
      plot = surv_plot$plot,
      width = 6, height = 5, dpi = 300
    )
    
    return(data.frame(
      gene = gene,
      pval = pval,
      survival_direction = survival_direction,
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}

# Run and save only significant results
results_LUAD <- do.call(rbind, lapply(sig_LUAD_genes, function(g) {
  run_survival_analysis(g, LUAD_data$expr, LUAD_data$clinical, "LUAD")
}))

results_LUSC <- do.call(rbind, lapply(sig_LUSC_genes, function(g) {
  run_survival_analysis(g, LUSC_data$expr, LUSC_data$clinical, "LUSC")
}))

# Save filtered results
write.csv(results_LUAD, "LUAD_significant_survival_summary.csv", row.names = FALSE)
write.csv(results_LUSC, "LUSC_significant_survival_summary.csv", row.names = FALSE)
