###======================== RESEARCH PROJECT====================================

#-------------------------------------------------------------------------------

#----------------DOWNLOADING DATA FOR LUAD AND LUSC-----------------------------
#Download tumor data for LUAD and LUSC that will contain the gene expression counts for each gene

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
  workflow.type = "STAR - Counts",  # STAR alignment workflow for count data
  sample.type = "Primary Tumor"  # Only include primary tumor samples
)

# Query for LUSC RNA-Seq data
query_LUSC <- GDCquery(
  project = "TCGA-LUSC",  
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"  
)

# Download the RNA-Seq data for LUAD and LUSC
GDCdownload(query_LUAD)
GDCdownload(query_LUSC)

# Prepare the downloaded data for analysis 
# converts raw GDC format to summarized experiment object containing:
# gene expression matrix, sample metadata, gene metadata
LUAD_data <- GDCprepare(query_LUAD)
LUSC_data <- GDCprepare(query_LUSC)

#---------------------------DOWNLOAD lncRNA Data -------------------------------
#Extracting only lncRNA data from the dataset downloaded from TCGA using Esnembl annotations

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
#Get all annotated human genes
all_genes <- genes(EnsDb.Hsapiens.v86)

# Check available gene biotypes to identify relevant categories for lncRNAs
unique(all_genes$gene_biotype)

# Define a list of biotypes corresponding to long non-coding RNAs (lncRNAs)
lncRNA_types <- c("lincRNA","macro_lncRNA","bidirectional_promoter_lncRNA",
                  "sense_intronic", "sense_overlapping", "antisense")


# Extract only the genes classified as lncRNAs
lncRNA_genes <- all_genes[all_genes$gene_biotype %in% lncRNA_types, ]

# Print a few lncRNA gene IDs to verify extraction
head(lncRNA_genes$gene_id)

# Ensure gene IDs in RNA-seq data match Ensembl format 
#Find the first dot in the gene ID and remove everything after it
rownames(LUAD_data) <- sub("\\..*", "", rownames(LUAD_data))
rownames(LUSC_data) <- sub("\\..*", "", rownames(LUSC_data))

#MATCHING lncRNAs IN DATA
# Find lncRNA genes that are present in LUAD and LUSC RNA-Seq data
lncRNA_ids_LUAD <- intersect(rownames(LUAD_data), lncRNA_genes$gene_id)
lncRNA_ids_LUSC <- intersect(rownames(LUSC_data), lncRNA_genes$gene_id)

# Check the number of matched lncRNAs in both datasets
# This will give you the count of lncRNA genes found in each dataset
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
#Downloading the same data for normal samples in order to normalize counts and compare tumor vs normal
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


#------------------------ Prepare datasets for normalization -------------------
#Load Required Libraries 
library(DESeq2)    #Analyzing count data
library(tibble)    #Handle dataframes
library(readr)     #Read CSV files

# Loading data
# Read in tumor and control lncRNA counts for LUAD and LUSC
LUAD_counts <- read.csv("LUAD_lncRNA.csv", row.names = 1)
LUAD_controls <- read.csv("LUAD_control_lncRNA.csv", row.names = 1)

LUSC_counts <- read.csv("LUSC_lncRNA.csv", row.names = 1)
LUSC_controls <- read.csv("LUSC_control_lncRNA.csv", row.names = 1)

# Fix sample names (TCGA usually uses '-' not '.') so samples will match
#gsub (): global substitution that replaces literal dots with hyphens
#e.g. TCGA.AB.1234-->TCGA-AB-1234
fix_colnames <- function(df) {
  colnames(df) <- gsub("\\.", "-", colnames(df))
  return(df)
}

LUAD_counts <- fix_colnames(LUAD_counts)
LUAD_controls <- fix_colnames(LUAD_controls)
LUSC_counts <- fix_colnames(LUSC_counts)
LUSC_controls <- fix_colnames(LUSC_controls)

# Add suffix to each sample to distinguish the tumor and normal samples
#paste0(): concnatinates strings --> adds "_tumor" or "_normal" to each column name
make_unique_colnames <- function(df, label) {
  colnames(df) <- paste0(colnames(df), "_", label)
  return(df)
}


LUAD_counts <- make_unique_colnames(LUAD_counts, "tumor")
LUAD_controls <- make_unique_colnames(LUAD_controls, "normal")
LUSC_counts <- make_unique_colnames(LUSC_counts, "tumor")
LUSC_controls <- make_unique_colnames(LUSC_controls, "normal")


#Get the number of samples
cat("LUAD Tumor samples:", ncol(LUAD_counts), "\n")
cat("LUAD Normal samples:", ncol(LUAD_controls), "\n")
cat("LUSC Tumor samples:", ncol(LUSC_counts), "\n")
cat("LUSC Normal samples:", ncol(LUSC_controls), "\n")

# Combine tumor and control counts 
LUAD_combined <- cbind(LUAD_counts, LUAD_controls)
LUSC_combined <- cbind(LUSC_counts, LUSC_controls)

# Save for future use if needed
write.csv(LUAD_combined, "LUAD_combined.csv")
write.csv(LUSC_combined, "LUSC_combined.csv") 

#=======Check if normal and tumor tissues come from the same paient

# Function to extract patient ID
get_patient_id <- function(barcode) {
  # Remove _tumor/_normal
  barcode <- sub("_(tumor|normal)$", "", barcode)
  # Replace dots with dashes
  barcode <- gsub("\\.", "-", barcode)
  # Split into parts
  parts <- strsplit(barcode, "-")[[1]]
  # Return first 3 parts
  paste(parts[1:3], collapse = "-")
}

# Apply to LUAD and LUSC
luad_tumor_ids  <- sapply(luad_tumor, get_patient_id)
luad_normal_ids <- sapply(luad_normal, get_patient_id)
lusc_tumor_ids  <- sapply(lusc_tumor, get_patient_id)
lusc_normal_ids <- sapply(lusc_normal, get_patient_id)

# Build match tables
luad_matches <- data.frame(
  Normal_Sample = luad_normal,
  Patient_ID = luad_normal_ids,
  Matched_To_Tumor = luad_normal_ids %in% luad_tumor_ids
)

lusc_matches <- data.frame(
  Normal_Sample = lusc_normal,
  Patient_ID = lusc_normal_ids,
  Matched_To_Tumor = lusc_normal_ids %in% lusc_tumor_ids
)

# Check
table(luad_matches$Matched_To_Tumor)
table(lusc_matches$Matched_To_Tumor)



#------------------------------- Normalize counts ------------------------------

#Create a codition label for DESeq2
#Create a factor that tells whether each sample is a tumor or normal sample
LUAD_condition <- factor(c(rep("Tumor", ncol(LUAD_counts)), rep("Normal", ncol(LUAD_controls))))
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



#=====Before and after VST plot and normalization for LUAD and LUSC

# Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)

#Function for processing and plotting
plot_expression <- function(raw_file, norm_file, vst_file, cancer_type) {
  
  # Load raw, normalized vst expression data
  raw <- read.csv(raw_file, row.names = 1)
  norm <- read.csv(norm_file, row.names = 1)
  vst <- read.csv(vst_file, row.names = 1)
  
  # Keep only tumor samples
  tumor_samples <- grep("_tumor$", colnames(raw), value = TRUE)
  
  # Select 30 random tumor samples
  # Set.seed to get the same result everytime we run the code
  set.seed(310)
  selected_samples <- sample(tumor_samples, 30)
  
  # Subset data using only the 30 selected samples
  raw_subset <- raw[, selected_samples]
  vst_subset <- vst[, selected_samples]
  
  # Melt to long format for plotting
  raw_long <- melt(raw_subset, variable.name = "Sample", value.name = "Counts")
  vst_long <- melt(vst_subset, variable.name = "Sample", value.name = "Counts")
  
  # Label stage
  raw_long$Stage <- "Raw"
  vst_long$Stage <- "VST"
  
  # Combine the data from raw and vst data
  combined <- rbind(raw_long, vst_long)
  
  # Add cancer type to the data
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
# This step stabilizes variance estimates for lowly expressed genes
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


#----------------------------MA plot to visualize DEA---------------------------
#Load library
library(tidyverse)

# Load DESeq2 results 
res_LUAD <- read.csv("LUAD_DESeq2_all_results.csv", row.names = 1)
res_LUSC <- read.csv("LUSC_DESeq2_all_results.csv", row.names = 1)

# Create an MA plot
# x-axis= baseMean from the DESeq2 results (mean expression)
# y-axis= log2foldchange from the DESeq2 results ( change in expression between tumor vs normal)
#color code points based on significance --> padj <0.05 significant=TRUE=deep red or else light red
# Use log10 scale in the x axis to spread out the low expression genes

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


#-------------------Venn to find shared and unique DEGS-------------------------
#Load Libraries 
library(dplyr)  #Data manipulation
library(VennDiagram)  #Draw a venn diagram


# Load DESeq2 shrunken results (already filtered for p-value & log2FC)
res_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
res_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)


# Extract gene IDS which will be the input for the Venn diagram
sig_LUAD <- rownames(res_LUAD)
sig_LUSC <- rownames(res_LUSC)

#Create the venn
venn <- venn.diagram(
  #Create a two-set venn digram comparing LUAD vs LUSC
  x = list(
    LUAD = sig_LUAD,
    LUSC = sig_LUSC
  ),
  #Title groups
  category.names = c("LUAD DEGs", "LUSC DEGs"),
  #Save the diagram
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
  
  # Font sizes and bold styling 
  cex = 4,                  # Count size
  cat.cex = 4,              # Category label size
  main.cex = 3.9,             # Title size
  fontface = "bold",        # Counts bold
  cat.fontface = "bold",    # Category labels bold
  main.fontface = "bold",   # Title bold
  #Use sans font --> clean font style
  fontfamily = "sans",         #Font family for numbers inside the circles
  cat.fontfamily = "sans",     #Font family for category labels ( LUAD DEGs and LUSC DEGs)
  main.fontfamily = "sans",    #Font family for the main title
  
  # Label positioning
  cat.pos = c(-155, 155),
  cat.dist = c(0.05, 0.05),
  
  # Title
  main = "SHARED AND UNIQUE DE lncRNAs"
)




#------------------- Map Ensembl IDs to Gene Names using biomaRt ---------------
#Load libraries
library(biomaRt)


# Load DESeq2 shrunken results (already filtered for p-value & log2FC)
deg_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
deg_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)


# Connect to Esnemble BioMart server to get the annotations
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract all unique Ensembl IDs from LUAD and LUSC
all_ensembl_ids <- unique(c(rownames(deg_LUAD), rownames(deg_LUSC)))


# Query gene names
# Use biomaRt to map Ensembl IDs to HGNC symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_ensembl_ids,
  mart = mart
)


# Save gene map
write.csv(gene_map, "gene_map_LUAD_LUSC.csv", row.names = FALSE)

# Create master gene label vector for quick renaming during plotting
master_gene_label <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)

# Print a sample of the mapping
print(head(master_gene_label))

# Save to RData file for future use
save(master_gene_label, file = "master_gene_label.RData")



#---------------------------- Select Top 15 DEGs -------------------------------
#Load libraries
library(dplyr)

# Load DEGs
deg_LUAD <- read.csv("LUAD_DESeq2_all_results_shrunk.csv", row.names = 1)
deg_LUSC <- read.csv("LUSC_DESeq2_all_results_shrunk.csv", row.names = 1)

# Filter significant DEGs first based on absLog2FC and adjusted p-value
# Arrange based on lowest padj and highest absLog2FC
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


#-------------------------------BOXPLOTS OF TOP 15 GENES -----------------------
# Load required libraries
library(dplyr) #Data manipulation
library(ggplot2)  #Plotting
library(reshape2)  #Reshaping the data
library(ggpubr)  #Statistical tests for ggplot2


# Load log-transformed data and preserve rownames as a column called 'Gene'
LUAD_vst <- read.csv("LUAD_lncRNA_log_trans_after_norm.csv") |> 
  rename(Gene = X)
LUSC_vst <- read.csv("LUSC_lncRNA_log_trans_after_norm.csv") |> 
  rename(Gene = X)


# Filter to Top 15 genes from the log transformed count
# Fix also Ensemble IDs by removing the version numbers
LUAD_top <- LUAD_vst[LUAD_vst$Gene %in% gsub("\\..*", "", top15_LUAD), ]
LUSC_top <- LUSC_vst[LUSC_vst$Gene %in% gsub("\\..*", "", top15_LUSC), ]

# Melt to long format
#Columns: Gene=Ensembl IDs, Variable=Sample IDs, Value=Expression levels
LUAD_melt <- melt(LUAD_top, id.vars = "Gene")
LUSC_melt <- melt(LUSC_top, id.vars = "Gene")


# Add SampleType column based on name suffix
LUAD_melt$SampleType <- ifelse(grepl("_tumor$", LUAD_melt$variable, ignore.case = TRUE), "Tumor", "Normal")
LUSC_melt$SampleType <- ifelse(grepl("_tumor$", LUSC_melt$variable, ignore.case = TRUE), "Tumor", "Normal")


# Replace Gene Ensembl ID with gene symbol from the map we created before
# Ensembl ID if missing then Ensembl ID will be used
LUAD_melt$Gene <- ifelse(
  #If the gene symbol is missing or empty, use the Ensembl ID
  is.na(master_gene_label[LUAD_melt$Gene]) | master_gene_label[LUAD_melt$Gene] == "",
  LUAD_melt$Gene,
  #If the gene symbol is available, use it
  master_gene_label[LUAD_melt$Gene]  
)

LUSC_melt$Gene <- ifelse(
  is.na(master_gene_label[LUSC_melt$Gene]) | master_gene_label[LUSC_melt$Gene] == "",
  LUSC_melt$Gene,
  master_gene_label[LUSC_melt$Gene]  
)

#Create boxplots for top15 DEGs comparing tumor vs normal VST
# x-axis = Gene names
# y-axis = VST expression
# fill= SampleType (Tumor vs Normal)
# stat_compare_means --> adds p-values from Wilcoxon test comparing normal vs tissue


# Plot LUAD
p_luad <- ggplot(LUAD_melt, aes(x = Gene, y = value, fill = SampleType)) +
  geom_boxplot(width = 0.85,outlier.shape = NA, alpha = 0.7) +  
  stat_compare_means(aes(group = SampleType),
                     method = "wilcox.test",label = "p.signif",size = 12, tip.length = 0.01, label.y = max(LUAD_melt$value) * 1.2) +
  scale_fill_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
  coord_cartesian(ylim = c(min(LUAD_melt$value), max(LUAD_melt$value) * 1.3)) +
  theme_minimal(base_size = 18) +
  ggtitle("Top 15 Differentially Expressed lncRNAs in LUAD") +
  ylab("VST") +xlab("Gene") +
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



ggsave("LUAD_Top15_Boxplot_with_GeneNames.png", plot = p_luad, width = 20, height = 10, dpi = 300)


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
ggsave("LUSC_Top15_Boxplot_with_GeneNames.png", plot = p_lusc, width = 17, height = 10, dpi = 300)

# Show both plots
print(p_luad)
print(p_lusc)




#----------------------- Download clinical data  -------------------------------
# Load Required Libraries
library(TCGAbiolinks) #Downloading and handling TCGA data
library(dplyr) #Data cleaning, filtering, selecting, ordering
library(readr) #Reading CSV files

# Load clinical metadata from TCGA
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD")
clinical_LUSC <- GDCquery_clinic(project = "TCGA-LUSC")

# Convert vital_status into binary outcome
# Create deceased column that indicates if the patient is deceased
clinical_LUAD$deceased <- ifelse(clinical_LUAD$vital_status == "Alive", FALSE, TRUE)
clinical_LUSC$deceased <- ifelse(clinical_LUSC$vital_status == "Alive", FALSE, TRUE)

# Create overall_survival column  that will indicate the survival time in days
# If deceased, use days_to_death; otherwise, use days_to_last_follow_up
clinical_LUAD$overall_survival <- ifelse(clinical_LUAD$deceased,
                                         clinical_LUAD$days_to_death,
                                         clinical_LUAD$days_to_last_follow_up)

clinical_LUSC$overall_survival <- ifelse(clinical_LUSC$deceased,
                                         clinical_LUSC$days_to_death,
                                         clinical_LUSC$days_to_last_follow_up)


# Save cleaned clinical data to new variables
clinical_LUAD_clean <- clinical_LUAD
clinical_LUSC_clean <- clinical_LUSC


# Write to CSV files if needed
write_csv(clinical_LUAD_clean, "clinical_LUAD_clean.csv")
write_csv(clinical_LUSC_clean, "clinical_LUSC_clean.csv")


#-----------------Summary for data: number of samples, gender etc---------------
# Load libraries
library(dplyr)   #For data manipulation
library(readr)   #Reading csvs
library(kableExtra)    #For creating a clean LaTex table


# Load data
expr_LUAD <- read.csv("LUAD_combined.csv", row.names = 1)
expr_LUSC <- read.csv("LUSC_combined.csv", row.names = 1)

clinical_LUAD <- read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)


# Helper to clean expression sample IDs
# This function will remove the "_tumor" or "_normal" suffix and replace "." with "-" in the sample IDs
clean_expr_ids <- function(expr) {
  ids <- gsub("_tumor$|_normal$", "", colnames(expr))
  ids <- gsub("\\.", "-", ids)
  sapply(strsplit(ids, "-"), function(x) paste(x[1:3], collapse = "-"))
}
#Apply function
expr_ids_LUAD <- clean_expr_ids(expr_LUAD)
expr_ids_LUSC <- clean_expr_ids(expr_LUSC)


#Extract sample IDs by keeping the first 12 characters
clinical_LUAD$short_id <- substr(clinical_LUAD$submitter_id, 1, 12)
clinical_LUSC$short_id <- substr(clinical_LUSC$submitter_id, 1, 12)

#Match samples between expression and clinical data
matched_LUAD <- clinical_LUAD[clinical_LUAD$short_id %in% expr_ids_LUAD, ]
matched_LUSC <- clinical_LUSC[clinical_LUSC$short_id %in% expr_ids_LUSC, ]

# Helper to clean clinical data
clean_clinical <- function(df) {
  #Convert age to years from days
  df$age_years <- df$age_at_diagnosis / 365.25
  # Extract the roman numeral from the AJCC stages
  df$stage_grouped <- gsub("Stage ([IV]+)[A-B]?", "\\1", df$ajcc_pathologic_stage)
  # Convert stage to factor with ordered levels
  df$stage_grouped <- factor(df$stage_grouped, levels = c("I", "II", "III", "IV"))
  #Group smoking status into categories
  df$smoking_group <- case_when(
    #Never smokers--> if it contains "non-smoker" or "never"
    grepl("(?i)non-smoker|never", df$tobacco_smoking_status) ~ "Never Smoker",
    # Current smokers --> if it contains "current smoker"
    grepl("(?i)current smoker", df$tobacco_smoking_status) ~ "Current Smoker",
    #Ex-smokers --> if it contains "ex-smoker" or "reformed"
    grepl("(?i)reformed", df$tobacco_smoking_status) ~ "Ex-Smoker",
    # If it contains "unknown" or is NA, label as "Unknown"
    TRUE ~ "Unknown"
  )
  return(df)
}


#Apply function 
matched_LUAD <- clean_clinical(matched_LUAD)
matched_LUSC <- clean_clinical(matched_LUSC)



# Get tumor and normal counts
get_sample_counts <- function(expr) {
  #Get column names that end with _tumor and get their sum
  tumor <- sum(grepl("_tumor$", colnames(expr)))
  #Get column names that end with _normal and get their sum
  normal <- sum(grepl("_normal$", colnames(expr)))
  #The result should be number of tumor/number of normal format
  paste0(tumor, " / ", normal)
}

#Apply function
tumor_normal_LUAD <- get_sample_counts(expr_LUAD)
tumor_normal_LUSC <- get_sample_counts(expr_LUSC)


# Build summary table
summary_table <- tibble(
  #Create a table with three columns: Variable, LUAD and LUSC
  #Variable column includes number of samples in expression data and samples that we have demographic information about
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
    #Counts of tumor/counts of normal
    tumor_normal_LUAD,
    #Number of expression sampls that also have clinical data
    nrow(matched_LUAD),
    #Age at diagnosis--> Mean and Standard deviation of age in years
    sprintf("%.1f ± %.1f", mean(matched_LUAD$age_years, na.rm = TRUE), sd(matched_LUAD$age_years, na.rm = TRUE)),
    #Minimum and maximum age at diagnosis
    sprintf("%.1f – %.1f", min(matched_LUAD$age_years, na.rm = TRUE), max(matched_LUAD$age_years, na.rm = TRUE)),
    #Counts how many males/females in the matched data
    paste(table(matched_LUAD$gender), collapse = " / "),
    #Counts how many current/ex smokers/never smokers there are in matched data
    paste(table(matched_LUAD$smoking_group)[c("Current Smoker", "Ex-Smoker", "Never Smoker")], collapse = " / "),
    #Counts how many patients there are in each cancer stage
    paste(table(matched_LUAD$stage_grouped), collapse = " / "),
    #Counts how mant people are alive/dead
    paste(table(matched_LUAD$vital_status), collapse = " / ")
  ),
  
  #Similarly for LUSC
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
#kbl: creates a table using the summary we created
#kable_Styling--> striped (striped (alternating row colors), hover(highlighting row on mouse hover), condensed (more compact spacing), full_width=FALSE (keeps table compact and centerd))
kbl(summary_table) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))




# -------- Build tumor-only, clinically matched expression matrices ------------
#Load library
library(readr)


# helper: TCGA submitter_id (TCGA-XX-XXXX) from column names like "TCGA-XX-XXXX_tumor"
to_submitter3 <- function(x) {
  x <- gsub("(_tumor|_normal)$", "", x, ignore.case = TRUE)  # drop suffix
  x <- gsub("\\.", "-", x)                                   # dots -> dashes
  vapply(strsplit(x, "-"), function(p) {
    if (length(p) >= 3) paste(p[1:3], collapse = "-") else NA_character_
  }, character(1))
}

# collapse duplicate columns (replicates) by mean
collapse_dupes_by_mean <- function(mat) {
  ids <- colnames(mat)
  idx <- split(seq_along(ids), ids)
  out <- lapply(idx, function(j)
    if (length(j) == 1) as.numeric(mat[, j])
    else rowMeans(as.matrix(mat[, j, drop = FALSE]), na.rm = TRUE)
  )
  out <- do.call(cbind, out)
  colnames(out) <- names(idx); rownames(out) <- rownames(mat)
  out
}

# generic: from any expression CSV -> tumor-only, clinically matched matrix
make_tumor_expr_matched <- function(expr_file, clinical_file, out_expr_csv, cancer_tag = "") {
  
  # load expression & clinical
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)  
  clin <- read_csv(clinical_file, show_col_types = FALSE)
  
  # keep tumor-only columns (end with "_tumor")
  tumor_cols <- grep("_tumor$", colnames(expr), value = TRUE, ignore.case = TRUE)
  expr_tumor <- expr[, tumor_cols, drop = FALSE]
  
  # harmonize colnames to submitter_id format (TCGA-XX-XXXX) and collapse duplicates
  colnames(expr_tumor) <- to_submitter3(colnames(expr_tumor))
  expr_tumor <- expr_tumor[, !is.na(colnames(expr_tumor)), drop = FALSE]
  expr_tumor <- collapse_dupes_by_mean(expr_tumor)
  
  # strip Ensembl version if present
  rownames(expr_tumor) <- gsub("\\..*", "", rownames(expr_tumor))
  
  # align to clinical submitter_id (first 12 chars)
  clin_ids <- substr(clin$submitter_id, 1, 12)
  keep_ids <- intersect(colnames(expr_tumor), clin_ids)
  
  expr_matched <- expr_tumor[, keep_ids, drop = FALSE]
  # (optional) reorder clinical if you’ll use it later side-by-side:
  # clin <- clin[match(keep_ids, clin_ids), ]
  
  # save
  write.csv(expr_matched, out_expr_csv)
  cat("[", cancer_tag, "] tumor-only matched samples:", length(keep_ids), "\n")
}

#Apply function
make_tumor_expr_matched(
  expr_file = "LUAD_lncRNA_normalized_counts_DESeq2.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  out_expr_csv = "LUAD_expr_matched.csv",
  cancer_tag = "LUAD"
)

make_tumor_expr_matched(
  expr_file = "LUSC_lncRNA_normalized_counts_DESeq2.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  out_expr_csv = "LUSC_expr_matched.csv",
  cancer_tag = "LUSC"
)






#-----------------------PCA for LUAD and LUSC lncRNAs---------------------------
# Load required libraries
library(ggplot2)
library(DESeq2)
library(dplyr)

# Function to run PCA and plot
run_pca_plot <- function(vst_file, cancer_type) {
  
  # Load vst data
  vst_data <- read.csv(vst_file, row.names = 1)
  
  # Transpose: genes as columns, samples as rows
  vst_t <- t(vst_data)
  
  # Detect if a sample is tumor or normal
  sample_type <- ifelse(grepl("_tumor$", rownames(vst_t)), "Tumor", "Normal")
  #Save this info into a table
  sample_df <- data.frame(Sample = rownames(vst_t), Type = sample_type)
  
  # Run PCA
  pca <- prcomp(vst_t, scale. = TRUE)
  # Create a data frame for PCA results
  pca_df <- as.data.frame(pca$x)
  # Add sample type to PCA results
  pca_df$Type <- sample_type
  
  # Calculate percentage of variance explained by each PC
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  #Plot PC1 vs PC2
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Type)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(
      title = paste("PCA -", cancer_type),
      x = paste0("PC1 (", percentVar[1], "% variance)"),
      y = paste0("PC2 (", percentVar[2], "% variance)")
    ) +
    scale_color_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
}

# Run and plot for LUAD
p1 <- run_pca_plot("LUAD_lncRNA_log_trans_after_norm.csv", "LUAD")
ggsave("PCA_LUAD.png", plot = p1, width = 8, height = 6, dpi = 300)
print(p1)

# Run and plot for LUSC
p2 <- run_pca_plot("LUSC_lncRNA_log_trans_after_norm.csv", "LUSC")
ggsave("PCA_LUSC.png", plot = p2, width = 8, height = 6, dpi = 300)
print(p2)



#----------------Early detection: Boxplots of Top 15 Genes by Stage ------------
#Load the necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)


# Function to plot stage vs normal for top 15 genes
plot_stage_vs_normal_top15_facet_singlefile <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  # Load expression data and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Split sample types based on tumor or normal suffix
  tumor_cols <- grep("_tumor$", colnames(expr), value = TRUE)
  normal_cols <- grep("_normal$", colnames(expr), value = TRUE)
  
  # Extract tumor and normal expression matrices
  expr_tumor <- expr[, tumor_cols]
  expr_normal <- expr[, normal_cols]
  
  # Fix sample IDs 
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(x, split = "\\."))
      parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  
  # Apply the function to fix sample IDs
  colnames(expr_tumor) <- fix_sample_ids(tumor_cols)
  colnames(expr_normal) <- fix_sample_ids(normal_cols)
  
  # Match tumor samples to clinical
  tumor_samples <- colnames(expr_tumor)
  clinical_match <- clinical[match(tumor_samples, clinical$submitter_id), ]
  
  #Extract stage information
  stage_group <- gsub("^Stage\\s*([IV]+)[A-Z]*$", "Stage \\1", clinical_match$ajcc_pathologic_stage)
  
  #Assign Stage I and II groups
  group_tumor <- ifelse(stage_group %in% c("Stage I", "Stage II"), stage_group, NA)
  
  # Assign Normal group
  group_normal <- rep("Normal", length(normal_cols))
  
  # Get top 15 genes
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Filter expression matrices for top 15 genes
  expr_tumor <- expr_tumor[rownames(expr_tumor) %in% top15_genes_clean, , drop = FALSE]
  expr_normal <- expr_normal[rownames(expr_normal) %in% top15_genes_clean, , drop = FALSE]
  
  # Combine genes from both tumor and normal
  genes_to_plot <- unique(c(rownames(expr_tumor), rownames(expr_normal)))
  
  # Reshape data to long format
  all_gene_data <- list()
  
  # Loop through each gene and create a data frame for plotting
  for (gene in genes_to_plot) {
    #Get gene symbol from master gene label
    gene_symbol <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    # Tumor
    if (gene %in% rownames(expr_tumor)) {
      # Extract expression vector for the gene in tumor samples
      expr_vec_tumor <- as.numeric(expr_tumor[gene, ])
      # Create a data frame for tumor samples with expression and group (Stage I/II) and gene
      df_tumor <- data.frame(
        expression = expr_vec_tumor,
        group = group_tumor,
        Gene = gene_symbol
      )
      # Add to the list of all gene data
      all_gene_data[[paste0(gene, "_tumor")]] <- df_tumor
    }
    
    # Normal
    if (gene %in% rownames(expr_normal)) {
      # Extract expression vector for the gene in normal samples
      expr_vec_normal <- as.numeric(expr_normal[gene, ])
      # Create a data frame for normal samples with expression and group (Normal) and geneSymbol
      df_normal <- data.frame(
        expression = expr_vec_normal,
        group = group_normal,
        Gene = gene_symbol
      )
      # Add to the list of all gene data
      all_gene_data[[paste0(gene, "_normal")]] <- df_normal
    }
  }
  
  # Combine all gene data frames into one
  df <- bind_rows(all_gene_data)
  # Remove rows with NA in group
  df <- df[!is.na(df$group), ]
  # Convert group to factor with specific levels
  df$group <- factor(df$group, levels = c("Normal", "Stage I", "Stage II"))
  
  # For each gene and group show how many samples contributed to the plot
  label_df <- df %>%
    group_by(Gene, group) %>%
    summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.12, .groups = "drop")
  
  # Define comparisons to be used in the wilcoxon test
  comparisons_list <- list(
    c("Normal", "Stage I"),
    c("Normal", "Stage II")
  )
  
  # Plot
  p <- ggplot(df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_text(data = label_df, aes(x = group, y = y, label = paste0("n = ", n)),
              inherit.aes = FALSE, size = 9, color = "black", fontface = "bold") +
    stat_compare_means(comparisons = comparisons_list, method = "wilcox.test",label = "p.signif",size = 12,tip.length = 0.015, step.increase = 0.15, bracket.size = 2,bracket.colour = "black") +
    scale_fill_manual(values = c("Normal" = "#F4A8A8", "Stage I" = "#A0A0A0", "Stage II" = "#A80000")) +
    facet_wrap(~Gene, ncol = 4, scales = "free_y") +
    coord_cartesian(clip = "off") +
    labs(title = paste0("Top 15 lncRNA Expression by Stage (", cancer_type, ")"),
         x = "Group",y = "VST") +
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
  
  #Save the plot
  ggsave(paste0("Top15_lncRNAs_by_stage_faceted_", cancer_type, "_EARLY_ONLY_VISIBLE.png"),plot = p,
         width = 32, height = 36, dpi = 300)
  
  print(p)
  cat("Saved combined faceted plot for", cancer_type, "\n")
}


# Apply function for LUAD
plot_stage_vs_normal_top15_facet_singlefile(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,
  cancer_type = "LUAD"
)

# Apply function for LUSC
plot_stage_vs_normal_top15_facet_singlefile(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",
  clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,
  cancer_type = "LUSC"
)




# ------------------- General Survival Analysis --------------------------------
# Performs KM survival analysis for each top lncRNA without stratifying
library(survival)
library(survminer)
library(ggplot2)  # ggsave()
library(readr)    # Reading CSVs

#Re-read the matched expression files WITHOUT mangling names
# (check.names = FALSE keeps the hyphens in TCGA barcodes; otherwise R turns
# them into dots and sample IDs no longer match clinical submitter_id.)
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1, check.names = FALSE)
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1, check.names = FALSE)

#Just in case: normalize any dots that might still be there
# Ensure sample IDs look like "TCGA-XX-XXXX" (not "TCGA.XX.XXXX")
colnames(LUAD_expr) <- gsub("\\.", "-", colnames(LUAD_expr))
colnames(LUSC_expr) <- gsub("\\.", "-", colnames(LUSC_expr))

#Load clinical and coerce IDs to 12-char submitter_id
clinical_LUAD <- read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)

# Clinical IDs as 12-char submitter_id (TCGA-XX-XXXX) to match expression cols
clinical_LUAD$submitter_id <- substr(as.character(clinical_LUAD$submitter_id), 1, 12)
clinical_LUSC$submitter_id <- substr(as.character(clinical_LUSC$submitter_id), 1, 12)


#Function: KM + Cox per gene (High vs Low by median
combined_survival_and_cox <- function(expr_matrix, clinical_df, top15_genes,
                                      master_gene_label, cancer_type,
                                      save_all_plots = TRUE) {
  #  harmonize IDs & match samples (VERY IMPORTANT
  #  expr_matrix is genes x samples, with sample IDs as TCGA-XX-XXXX
  #  clinical_df has submitter_id; we match on the first 12 chars
  colnames(expr_matrix) <- gsub("\\.", "-", colnames(expr_matrix))  # safety
  clin_ids <- substr(clinical_df$submitter_id, 1, 12)
  
  # use only shared samples and put them in the SAME ORDER in both objects
  use_ids <- intersect(colnames(expr_matrix), clin_ids)
  if (length(use_ids) < 10) {
    stop(paste("Too few matched samples for", cancer_type, ":", length(use_ids)))
  }
  
  expr_matrix  <- expr_matrix[, use_ids, drop = FALSE]
  clinical_df  <- clinical_df[match(use_ids, clin_ids), , drop = FALSE]
  
  # quick sanity check
  if (ncol(expr_matrix) != nrow(clinical_df)) {
    stop("Internal matching error: ncol(expr_matrix) != nrow(clinical_df)")
  }
  
  # clean gene IDs (drop Ensembl version numbers like '.1', '.2', …)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  results <- list()
  
  # loop genes
  for (gene in top15_genes_clean) {
    # Check gene exists in the expression matrix
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    # Map Ensembl -> HGNC symbol (fallback to Ensembl ID if missing)
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene, master_gene_label[gene]
    )
    
    # Working clinical DF for this gene (ordering preserved)
    df <- clinical_df
    
    # Attach gene expression vector (length == nrow(df) by construction)
    df$expression <- as.numeric(expr_matrix[gene, use_ids])
    
    # Require survival columns to be present
    if (!("overall_survival" %in% names(df)) || !("deceased" %in% names(df))) {
      message("Missing survival columns for gene: ", gene_name)
      next
    }
    
    # Drop rows with missing expr/survival
    df <- df[!is.na(df$expression) & !is.na(df$overall_survival) & !is.na(df$deceased),
             , drop = FALSE]
    if (nrow(df) < 10) next
    
    # High/Low split using the median within the current gene-specific data
    df$group <- factor(
      ifelse(df$expression >= median(df$expression, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")  # Low as reference in Cox
    )
    if (length(unique(df$group)) < 2) next
    
    # Kaplan–Meier + log-rank
    fit  <- survfit(Surv(overall_survival, deceased) ~ group, data = df)
    test <- survdiff(Surv(overall_survival, deceased) ~ group, data = df)
    pval <- 1 - pchisq(test$chisq, df = 1)
    
    # Cox proportional hazards (High vs Low)
    cox <- tryCatch(
      summary(coxph(Surv(overall_survival, deceased) ~ group, data = df)),
      error = function(e) NULL
    )
    
    # Collect results if Cox ran
    if (!is.null(cox)) {
      results[[length(results) + 1]] <- data.frame(
        gene_id     = gene,
        gene_symbol = gene_name,
        HR          = round(cox$coefficients[, "exp(coef)"], 3),
        lower_CI    = round(cox$conf.int[, "lower .95"], 3),
        upper_CI    = round(cox$conf.int[, "upper .95"], 3),
        KM_pval     = signif(pval, 3),
        Cox_pval    = signif(cox$coefficients[, "Pr(>|z|)"], 3),
        n           = nrow(df)
      )
    }
    
    # Save plot if requested or significant by KM
    if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
      plot <- ggsurvplot(
        fit, data = df, pval = TRUE,
        title = paste(gene_name, "in", cancer_type),
        xlab = "Days", ylab = "Survival Probability"
      )
      ggsave(
        paste0(cancer_type, "_", gene_name, "_survival.png"),
        plot = plot$plot, width = 6, height = 5, dpi = 300
      )
      cat("Saved plot for gene:", gene_name, "\n")
    }
  }
  
  # Write results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_survival_results_ALL.csv"), row.names = FALSE)
    cat("Saved unstratified survival results for", nrow(results_df), "genes in", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}

# Run for LUAD and LUSC
combined_survival_and_cox(LUAD_expr, clinical_LUAD, top15_LUAD, master_gene_label, "LUAD", TRUE)
combined_survival_and_cox(LUSC_expr, clinical_LUSC, top15_LUSC, master_gene_label, "LUSC", TRUE)


#====== Plot Cox regression results (Forest Plot) ======

# Load required libraries
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation (e.g., filter, mutate, bind_rows)
library(readr)     # For reading CSV files


# Read the combined results for LUAD and LUSC from CSV files
luad <- read_csv("LUAD_survival_results_ALL.csv")
lusc <- read_csv("LUSC_survival_results_ALL.csv")


# Add a cancer type label and convert gene names to uppercase (for consistency)
luad <- luad %>% mutate(cancer = "LUAD", gene_symbol = toupper(gene_symbol))
lusc <- lusc %>% mutate(cancer = "LUSC", gene_symbol = toupper(gene_symbol))

# Combine LUAD and LUSC results into a single dataframe
combined <- bind_rows(luad, lusc)

# Filter genes of interest
# Select only specific genes and their cancer types (e.g., ZFPM2-AS1 in LUAD and LINC00511 in LUSC)
filtered <- combined %>%
  filter(
    (gene_symbol == "ZFPM2-AS1" & cancer == "LUAD") |
      (gene_symbol == "LINC00511" & cancer == "LUSC")
  ) %>%
  mutate(
    # Create a label combining gene name and cancer type for x-axis
    label = paste0(gene_symbol, " (", cancer, ")"),
    
    # Create a factor variable indicating statistical significance of Cox p-value
    significance = factor(
      ifelse(Cox_pval < 0.05, "Significant", "Not significant"),
      levels = c("Significant", "Not significant")  # Ensure consistent color mapping
    ),
    
    # Reverse label order for top-down plotting
    label = factor(label, levels = rev(label))
  )


# Create dummy legend data (for consistency)
# This creates invisible points that force ggplot to always show both legend categories,
# even if one of them is not present in the filtered data.
legend_df <- data.frame(
  significance = factor(c("Significant", "Not significant"),
                        levels = c("Significant", "Not significant")),
  x = Inf,  # Off-plot coordinates to make them invisible
  y = Inf
)

# Create the plot

ggplot(filtered, aes(x = label, y = HR, ymin = lower_CI, ymax = upper_CI)) +
  # Draw point-range for each gene: point = HR, range = 95% CI
  geom_pointrange(aes(color = significance), size = 1.2) +
  
  # Add dummy points for legend categories
  geom_point(
    data = legend_df,
    mapping = aes(x = x, y = y, color = significance),
    inherit.aes = FALSE,
    size = 4,
    show.legend = TRUE
  ) +
  
  # Add a horizontal line at HR = 1 (null hypothesis of no effect)
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  
  # Set custom colors for significance levels
  scale_color_manual(
    values = c("Significant" = "firebrick", "Not significant" = "grey60"),
    drop = FALSE  # Ensure both levels are shown even if absent
  ) +
  
  # Flip coordinates so genes are listed top-to-bottom on the y-axis
  coord_flip() +
  
  # Axis and legend labels
  labs(
    title = "Forest Plot of Selected Prognostic lncRNAs",
    x = "lncRNA (Cancer Type)",
    y = "Hazard Ratio (95% CI)",
    color = "Significance"
  ) +
  
  # Apply a clean theme with larger base font
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",                         # Place legend on top
    plot.title = element_text(face = "bold", hjust = 0.5)  # Center and bold title
  )



#------------------------- Correlate Age with Top 15 DEGs ----------------------
# Load required libraries 
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)


#Perform age-expression correlation for the 15 genes per cancer type
correlate_age <- function(expr_file, clinical_file, cancer_type, top15_genes, master_gene_label, sample_type = c("tumor", "normal")) {
  
  # Ensure sample_type is one of the allowed values
  sample_type <- match.arg(sample_type)  
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Filter expression columns based on sample type
  if (sample_type == "tumor") {
    expr <- expr[, grepl("_tumor$", colnames(expr))]
  } else {
    expr <- expr[, grepl("_normal$", colnames(expr))]
  }
  
  # Fix sample IDs to match clinical$submitter_id
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(x, split = "\\."))
      parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match clinical to expression samples
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Convert age from days to years
  clinical$age <- as.numeric(clinical$age_at_diagnosis) / 365
  
  # Run correlation for Top 15 DEGs
  results <- lapply(top15_genes, function(gene) {
    #If the gene is in the expression matrix, extract its expression values and age values
    if (gene %in% rownames(expr)) {
      expr_vals <- as.numeric(expr[gene, ])
      age_vals <- clinical$age
      
      # Check if both vectors have the same length
      if (length(expr_vals) == length(age_vals)) {
        # Perform Spearman correlation test
        test <- cor.test(expr_vals, age_vals, method = "spearman", exact = FALSE)
        # Return a data frame with the results that include gene ID, Spearman rho, and p-value
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
  
  #Check if results_df is NULL or empty
  if (is.null(results_df)) {
    cat("No valid correlations computed for", cancer_type, "-", sample_type, "\n")
    return(NULL)
  }
  
  # Adjust p-values using Benjamini-Hochberg method
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Clean gene IDs 
  results_df$gene_clean <- gsub("\\..*", "", results_df$gene_id)
  
  #Change gene IDs to gene symbols if they exist
  results_df$gene_symbol <- ifelse(
    is.na(master_gene_label[results_df$gene_clean]) | master_gene_label[results_df$gene_clean] == "",
    results_df$gene_clean,
    master_gene_label[results_df$gene_clean]
  )
  
  # Save results
  filename <- paste0(cancer_type, "_", sample_type, "_age_correlation_top15.csv")
  write.csv(results_df, filename, row.names = FALSE)
  cat("Done with", cancer_type, "-", sample_type, "- saved age correlation for Top 15 DEGs\n")
}

# Correlate age with LUAD samples
correlate_age("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, sample_type = "tumor")
correlate_age("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, sample_type = "normal")

# Correlate age with LUSC samples
correlate_age("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, sample_type = "tumor")
correlate_age("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, sample_type = "normal")


#===Plotting age correlation
#Load Required Libraries
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(grid)

#Load Correlation Results
# Load Correlation Results 
luad <- read.csv("LUAD_tumor_age_correlation_top15.csv")
lusc <- read.csv("LUSC_tumor_age_correlation_top15.csv")


#Format LUAD Matrix 
luad_filtered <- luad %>%
  #Filter out rows with NA in spearman_rho
  filter(!is.na(spearman_rho)) %>%
  # Mutate labels for significance
  #If padk<0.005 and abs spearman rho>0.1 then add a star to indicate a meaningful cor
  mutate(label = ifelse(padj < 0.05 & abs(spearman_rho) > 0.1,
                        paste0(sprintf("%.2f", spearman_rho), "*"),
                        sprintf("%.2f", spearman_rho))) %>%
  #Add gene symbols
  column_to_rownames("gene_symbol")


#Format LUSC Matrix as LUAD
lusc_filtered <- lusc %>%
  #Filter out rows with NA in spearman_rho
  filter(!is.na(spearman_rho)) %>%
  # Mutate labels for significance
  #If padk<0.005 and abs spearman rho>0.1 then add a star to indicate a meaningful cor
  mutate(label = ifelse(padj < 0.05 & abs(spearman_rho) > 0.1,
                        paste0(sprintf("%.2f", spearman_rho), "*"),
                        sprintf("%.2f", spearman_rho))) %>%
  #Add gene symbols
  column_to_rownames("gene_symbol")

#Extract the numerix spearman rho into a matric for plotting
luad_mat <- as.matrix(luad_filtered$spearman_rho)
lusc_mat <- as.matrix(lusc_filtered$spearman_rho)

#Ensure that each row of the correlation matrix has a corresponding gene symbol
rownames(luad_mat) <- rownames(luad_filtered)
rownames(lusc_mat) <- rownames(lusc_filtered)

#Creat label for the heatmap
luad_labels <- matrix(luad_filtered$label, ncol = 1)
lusc_labels <- matrix(lusc_filtered$label, ncol = 1)

#Assign gene symbols as row names to match the expression matrices
rownames(luad_labels) <- rownames(luad_filtered)
rownames(lusc_labels) <- rownames(lusc_filtered)

#Set column names for the matrices
colnames(luad_mat) <- "spearman_rho"
colnames(lusc_mat) <- "spearman_rho"

colnames(luad_labels) <- "spearman_rho"
colnames(lusc_labels) <- "spearman_rho"

# Create Side-by-Side Heatmaps with Larger Fonts
p1 <- pheatmap(luad_mat,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               display_numbers = luad_labels,   #Show correlation values
               main = "LUAD: Spearman Correlation with Age",
               color = colorRampPalette(c("#A80000", "white", "red"))(100),
               fontsize_row = 14,       # Increase gene name size
               fontsize_number = 12,    # Increase rho + star size
               silent = TRUE)

p2 <- pheatmap(lusc_mat,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               display_numbers = lusc_labels,
               main = "LUSC: Spearman Correlation with Age",
               color = colorRampPalette(c("#A80000", "white", "red"))(100),
               fontsize_row = 14,
               fontsize_number = 12,
               silent = TRUE)

# Save Combined Plot
png("Age_Correlation_Heatmap_LUAD_LUSC_with_stars.png", width = 2500, height = 1500, res = 200)  # also increase output size
grid.arrange(p1[[4]], p2[[4]], ncol = 2)
dev.off()

#==========Plot exp vs age bins
# Load required libraries
library(ggplot2)   # For plotting
library(ggpubr)    # For statistical comparisons (e.g., p-value annotations)
library(tidyverse) # For data wrangling

# Main function to plot expression of top genes across age bins, stratified by sample type (tumor/normal)
plot_age_bin_violin_per_gene_with_significance_flex_v5 <- function(vst_file, clinical_file, cancer_type, top15_genes, master_gene_label, pairwise = FALSE, sample_type = c("tumor", "normal")) {
  
  # Ensure sample_type is either "tumor" or "normal"
  sample_type <- match.arg(sample_type)
  
  # Load VST-normalized expression data and clinical metadata
  expr <- read.csv(vst_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Keep only columns (samples) matching the sample_type
  if (sample_type == "tumor") {
    expr <- expr[, grepl("_tumor$", colnames(expr))]
  } else {
    expr <- expr[, grepl("_normal$", colnames(expr))]
  }
  
  # Remove version suffix from Ensembl gene IDs
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Function to fix sample IDs to match clinical$submitter_id
  ffix_sample_ids <- function(ids) {
    ids <- gsub("(_tumor|_normal)$", "", ids, ignore.case = TRUE)  # remove suffix
    ids <- gsub("\\.", "-", ids)                                   # safety: dots -> dashes
    vapply(ids, function(x) {
      parts <- unlist(strsplit(x, "-"))
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA_character_
    }, character(1))
  }
  
  
  # Fix column names in expression matrix
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match clinical rows to expression columns
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Convert age from days to years
  clinical$age <- as.numeric(clinical$age_at_diagnosis) / 365
  
  # Loop through each gene and create a violin plot
  for (gene in top15_genes_clean) {
    
    # Skip if gene not found in expression data
    if (gene %in% rownames(expr)) {
      
      # Extract expression values for the gene
      expr_vals <- as.numeric(expr[gene, ])
      names(expr_vals) <- NULL  # Prevent warnings from unnamed numeric vectors
      
      # Create a dataframe combining expression and age
      df <- data.frame(
        expression = expr_vals,
        age = clinical$age
      )
      
      # Create custom age bins: 50–59, 60–69, 70–89 years
      df$age_bin <- cut(
        df$age,
        breaks = c(50, 60, 70, 90),
        right = FALSE,
        labels = c("50-59", "60-69", "70-89")
      )
      
      # Remove missing values and drop unused factor levels
      df <- df %>% dplyr::filter(!is.na(age_bin))
      df$age_bin <- droplevels(df$age_bin)
      
      # Skip gene if no valid samples remain
      if (nrow(df) == 0) {
        message("No valid data for gene ", gene, " — skipping.")
        next
      }
      
      # Get present age bins in data (may be <3 for some genes)
      bins_present <- sort(unique(df$age_bin))
      
      # Prepare pairwise comparisons if more than one bin is present
      my_comparisons <- if (length(bins_present) >= 2) {
        combn(as.character(bins_present), 2, simplify = FALSE)
      } else {
        NULL
      }
      
      # Convert Ensembl ID to gene symbol using master gene label (fallback to Ensembl ID if missing)
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Start building the plot
      p <- ggplot(df, aes(x = age_bin, y = expression)) +
        geom_violin(trim = FALSE, fill = "#A80000", color = "black", alpha = 0.7) +  # Violin plot
        geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5)                 # Boxplot overlay
      
      # Add statistical comparison (Wilcoxon pairwise or Kruskal-Wallis global)
      if (pairwise && !is.null(my_comparisons)) {
        p <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
          labs(subtitle = "Pairwise Wilcoxon tests")
      } else {
        p <- p + stat_compare_means(method = "kruskal.test", label = "p.signif") +
          labs(subtitle = "Global Kruskal-Wallis test")
      }
      
      # Finalize the plot with theme and axis labels
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
      
      # Save plot to PNG with filename indicating sample type and test type
      filename <- paste0(cancer_type, "_", gene_symbol, "_AgeBin_", toupper(sample_type), if (pairwise) "_PAIRWISE.png" else "_GLOBAL.png")
      ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
      cat("Saved plot:", filename, "\n")
      
    } else {
      message("Gene not found in expression matrix: ", gene)
    }
  }
  
  # Done message
  cat("Finished saving violin plots for", cancer_type, "-", sample_type, "\n")
}


# LUAD Tumor and Normal
plot_age_bin_violin_per_gene_with_significance_flex_v5("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, pairwise = TRUE, sample_type = "tumor")
plot_age_bin_violin_per_gene_with_significance_flex_v5("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, pairwise = TRUE, sample_type = "normal")

# LUSC Tumor and Normal
plot_age_bin_violin_per_gene_with_significance_flex_v5("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, pairwise = TRUE, sample_type = "tumor")
plot_age_bin_violin_per_gene_with_significance_flex_v5("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, pairwise = TRUE, sample_type = "normal")


#------------------Survival Analysis for DE lncRNAs correlated with age---------
# Load required libraries
library(survival)    # For survival models 
library(survminer)   # For visualization of survival curves with ggplot2

#Load required data
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1, check.names = FALSE)
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1, check.names = FALSE)
clinical_LUAD <- read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)

# This function performs survival analysis for each age group (<65 and ≥65 years)
# for each lncRNA significantly correlated with age.
# It saves Kaplan-Meier plots and returns only those with log-rank p < 0.05.
age_stratified_survival_and_cox <- function(expr_matrix, clinical_df, top15_genes,
                                            master_gene_label, cancer_type,
                                            save_all_plots = TRUE) {
  # Initialize results list
  results <- list()
  
  #ID HARMONIZATION & MATCHING (VERY IMPORTANT)
  # expr_matrix is genes x samples with TCGA-XX-XXXX as colnames.
  # clinical_df has submitter_id; we match on the first 12 chars.
  colnames(expr_matrix) <- gsub("\\.", "-", colnames(expr_matrix))  # safety
  clinical_df$submitter_id <- substr(as.character(clinical_df$submitter_id), 1, 12)
  
  use_ids <- intersect(colnames(expr_matrix), clinical_df$submitter_id)
  if (length(use_ids) < 10) {
    stop(paste("Too few matched samples for", cancer_type, ":", length(use_ids)))
  }
  
  # Keep only shared samples in the SAME ORDER in both objects
  expr_matrix <- expr_matrix[, use_ids, drop = FALSE]
  clinical_df <- clinical_df[match(use_ids, clinical_df$submitter_id), , drop = FALSE]
  
  # Quick sanity check
  if (ncol(expr_matrix) != nrow(clinical_df)) {
    stop("Internal matching error: ncol(expr_matrix) != nrow(clinical_df)")
  }
  
  #Clean gene IDs
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  # Loop through each gene
  for (gene in top15_genes_clean) {
    
    # Skip gene if not in the expression matrix
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    # Get gene symbol or fallback to Ensembl ID
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    # Build merged dataframe for this gene
    # Create merged dataframe
    df <- clinical_df
    # Attach expression vector (length equals nrow(df) by construction)
    df$expression <- as.numeric(expr_matrix[gene, use_ids])
    
    # Check required survival columns exist
    if (!all(c("overall_survival", "deceased", "age_at_diagnosis") %in% colnames(df))) {
      message("Missing required columns for gene: ", gene_name)
      next
    }
    
    # Convert age to years
    df$age <- as.numeric(df$age_at_diagnosis) / 365
    
    # Filter out missing values
    df <- df[!is.na(df$expression) & !is.na(df$overall_survival) &
               !is.na(df$deceased) & !is.na(df$age), , drop = FALSE]
    if (nrow(df) < 10) next
    
    # Assign age group
    df$age_group <- ifelse(df$age < 65, "<65", "≥65")
    
    # Loop over age groups
    for (age_grp in unique(df$age_group)) {
      
      # Subset data for this age group
      sub_df <- df[df$age_group == age_grp, , drop = FALSE]
      
      # Skip if too few samples
      if (nrow(sub_df) < 10) next
      
      # Split into High/Low expression based on median
      sub_df$group <- factor(
        ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High")   # Low = reference for Cox
      )
      
      # Skip if only one group present
      if (length(unique(sub_df$group)) < 2) next
      
      # Fit Kaplan-Meier curve
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      
      # Log-rank test
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Cox regression
      cox <- tryCatch({
        summary(coxph(Surv(overall_survival, deceased) ~ group, data = sub_df))
      }, error = function(e) NULL)
      
      # If Cox regression ran successfully, store results
      if (!is.null(cox)) {
        results[[length(results) + 1]] <- data.frame(
          gene_id   = gene,
          gene_symbol = gene_name,
          age_group = age_grp,
          HR        = round(cox$coefficients[, "exp(coef)"], 3),
          lower_CI  = round(cox$conf.int[, "lower .95"], 3),
          upper_CI  = round(cox$conf.int[, "upper .95"], 3),
          KM_pval   = signif(pval, 3),
          Cox_pval  = signif(cox$coefficients[, "Pr(>|z|)"], 3),
          n         = nrow(sub_df)
        )
      }
      
      # Save plot if required
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        # path-safe tag for filename
        age_tag <- ifelse(age_grp == "<65", "lt65", "ge65")
        
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = paste(gene_name, "in", cancer_type, "- Age", age_grp),
          xlab = "Days",
          ylab = "Survival Probability"
        )
        
        filename <- paste0(cancer_type, "_", gene_name, "_Age_", age_tag, "_survival.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        cat("Saved survival plot for:", gene_name, "- Age", age_grp, "\n")
      }
    }
  }
  
  # Combine and export results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_age_stratified_survival_results_ALL.csv"), row.names = FALSE)
    cat("Saved age-stratified survival results for", nrow(results_df),
        "gene/age group combinations in", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}


# Apply to LUAD
age_stratified_survival_and_cox(expr_matrix = LUAD_expr,clinical_df = clinical_LUAD,top15_genes = top15_LUAD,
                                master_gene_label = master_gene_label,cancer_type = "LUAD",save_all_plots = TRUE)

# Apply to LUSC
age_stratified_survival_and_cox(expr_matrix = LUSC_expr,clinical_df = clinical_LUSC,top15_genes = top15_LUSC,
                                master_gene_label = master_gene_label,cancer_type = "LUSC",save_all_plots = TRUE)


# --------------------- Correlate Gender with Top 15 DEGs ----------------------
# Load required libraries
library(dplyr)
library(readr)
library(rstatix)

# Perform gender-expression correlation for the Top 15 genes per cancer type
gender_correlation_top15 <- function(expr_file, clinical_file, cancer_type, top15_genes, master_gene_label, sample_type = c("tumor", "normal")) {
  
  # Ensure sample_type is one of the allowed values
  sample_type <- match.arg(sample_type)  
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Filter expression columns based on sample type
  if (sample_type == "tumor") {
    expr <- expr[, grepl("_tumor$", colnames(expr))]
  } else {
    expr <- expr[, grepl("_normal$", colnames(expr))]
  }
  
  # Fix sample IDs from the expression dataset to match submitter_id from the clinical dataset
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(x, split = "\\."))
      parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match clinical to expression samples
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Clean gene IDs in both the expression matrix and Top 15 list
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Run Wilcoxon rank-sum test for each Top 15 gene
  results <- lapply(top15_genes_clean, function(gene) {
    # If the gene is in the expression matrix, extract expression and gender values
    if (gene %in% rownames(expr)) {
      expr_vals <- as.numeric(expr[gene, ])
      gender_vals <- clinical$gender
      
      # Check that there are exactly 2 gender groups
      if (length(unique(na.omit(gender_vals))) == 2) {
        # Perform Wilcoxon rank-sum test
        test <- rstatix::wilcox_test(data.frame(expression = expr_vals, gender = gender_vals), expression ~ gender)
        # Calculate effect size
        effect <- rstatix::wilcox_effsize(data.frame(expression = expr_vals, gender = gender_vals), expression ~ gender)
        
        # Determine effect direction based on medians
        median_male <- median(expr_vals[gender_vals == "male"], na.rm = TRUE)
        median_female <- median(expr_vals[gender_vals == "female"], na.rm = TRUE)
        effect_direction <- ifelse(median_male > median_female, "male>female", "female>male")
        
        # Return results for this gene
        return(data.frame(
          gene_id = gene,
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
  
  # Combine results into a single data frame
  results_df <- do.call(rbind, results)
  
  # Check if results_df is NULL or empty
  if (is.null(results_df) || nrow(results_df) == 0) {
    cat("No valid gender comparisons for", cancer_type, "-", sample_type, "\n")
    return(NULL)
  }
  
  # Adjust p-values using Benjamini-Hochberg method
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Clean gene IDs for annotation
  results_df$gene_clean <- gsub("\\..*", "", results_df$gene_id)
  
  # Change gene IDs to gene symbols if available
  results_df$gene_symbol <- ifelse(
    is.na(master_gene_label[results_df$gene_clean]) | master_gene_label[results_df$gene_clean] == "",
    results_df$gene_clean,
    master_gene_label[results_df$gene_clean]
  )
  
  # Save results
  filename <- paste0(cancer_type, "_", sample_type, "_gender_correlation_top15.csv")
  write.csv(results_df, filename, row.names = FALSE)
  cat("Done with", cancer_type, "-", sample_type, "- saved gender correlation for Top 15 DEGs\n")
}

#Applying the function
# LUAD Tumor
gender_correlation_top15("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, sample_type = "tumor")

# LUAD Normal
gender_correlation_top15("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, sample_type = "normal")

# LUSC Tumor
gender_correlation_top15("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, sample_type = "tumor")

# LUSC Normal
gender_correlation_top15("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, sample_type = "normal")


#====== Plot for Gender-Correlated lncRNAs ======
# Load required libraries
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)


#Function that plots the expression distributions of each gender for each gene
plot_top15_gender_boxplots_for_given_genes <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type, sample_type = c("tumor", "normal")) {
  
  #Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Filter by sample type
  suffix <- ifelse(sample_type == "tumor", "_tumor", "_normal")
  expr <- expr[, grepl(suffix, colnames(expr))]
  
  # Fix sample IDs
  colnames(expr) <- gsub("(_tumor|_normal)$", "", colnames(expr))  # remove suffix
  colnames(expr) <- sapply(colnames(expr), function(x) {
    parts <- unlist(strsplit(x, split = "\\.")); parts <- parts[parts != ""]
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
  
  # Match clinical
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Clean gene IDs
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  #For each gene from the top15
  for (gene in top15_genes_clean) {
    
    #If it exists in the expression file
    if (!(gene %in% rownames(expr))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    #Map its ensembl id to gene symbol if it exists
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    #Create a dataframe that contains the expression values and the gender of each sample
    expr_vec <- as.numeric(expr[gene, ])
    df <- data.frame(expression = expr_vec, gender = clinical$gender)
    
    # Calculate the sample sizes
    label_df <- df %>%
      group_by(gender) %>%
      summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.10)
    
    #Create a list that states the comparison
    comparisons_list <- list(c("male", "female"))
    
    #Generate a plot that perfomrs wilcoxon test and states the significance by stars
    p <- ggplot(df, aes(x = gender, y = expression, fill = gender)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
      stat_compare_means(comparisons = comparisons_list, method = "wilcox.test",label = "p.signif", bracket.size = 0.5,step.increase = 0.1,size = 5) +
      geom_text(data = label_df,aes(x = gender, y = y, label = paste0("n = ", n)),
                inherit.aes = FALSE,size = 5.5,color = "black") +
      scale_fill_manual(values = c("male" = "grey", "female" = "#F8766D")) +
      labs(title = paste(gene_name, "expression in", cancer_type),
           x = "Gender", y = "Variance Stabilized Expression (VST)") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold")
      )
    
    ggsave(filename = paste0(cancer_type, "_", gene_name, "_gender_boxplot.png"),
           plot = p, width = 8, height = 6, dpi = 300
    )
    
    cat("Saved plot for gene:", gene_name, "\n")
  }
  
  cat("Finished plotting Top 15 gender boxplots for", cancer_type, "\n")
}



plot_top15_gender_boxplots_for_given_genes(expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUAD_clean.csv",top15_genes = top15_LUAD,
                                           master_gene_label = master_gene_label, cancer_type = "LUAD_Tumor", sample_type = "tumor")



plot_top15_gender_boxplots_for_given_genes(expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUAD_clean.csv",top15_genes = top15_LUAD,
                                           master_gene_label = master_gene_label,cancer_type = "LUAD_Normal",sample_type = "normal")


plot_top15_gender_boxplots_for_given_genes(expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUSC_clean.csv",top15_genes = top15_LUSC,
                                           master_gene_label = master_gene_label,cancer_type = "LUSC_Tumor",sample_type = "tumor")


plot_top15_gender_boxplots_for_given_genes(expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUSC_clean.csv",top15_genes = top15_LUSC,
                                           master_gene_label = master_gene_label,cancer_type = "LUSC_Normal",sample_type = "normal")




#--------------------------- Gender-Stratified Survival Analysis ---------------
# Load required libraries
library(survival)    # For survival models 
library(survminer)   # For visualization of survival curves with ggplot2


#Load required data
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1, check.names = FALSE)
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1, check.names = FALSE)
clinical_LUAD <- read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)


# This function performs survival analysis for each gender (male / female)
# for each lncRNA significantly correlated with gender.
# It saves Kaplan-Meier plots and returns results including KM p-value and Cox HR.
gender_stratified_survival_and_cox <- function(expr_matrix, clinical_df, top15_genes,
                                               master_gene_label, cancer_type,
                                               save_all_plots = TRUE) {
  
  # ID HARMONIZATION & MATCHING
  # expr_matrix is genes x samples with TCGA-XX-XXXX as colnames.
  # clinical_df has submitter_id; we match on the first 12 chars and REORDER
  # the clinical rows to be in the SAME order as the expression columns.
  colnames(expr_matrix) <- gsub("\\.", "-", colnames(expr_matrix))             # Just in case
  clinical_df$submitter_id <- substr(as.character(clinical_df$submitter_id), 1, 12)
  
  use_ids <- intersect(colnames(expr_matrix), clinical_df$submitter_id)
  if (length(use_ids) < 10) {
    stop(paste("Too few matched samples for", cancer_type, ":", length(use_ids)))
  }
  
  expr_matrix <- expr_matrix[, use_ids, drop = FALSE]
  clinical_df <- clinical_df[match(use_ids, clinical_df$submitter_id), , drop = FALSE]
  
  if (ncol(expr_matrix) != nrow(clinical_df)) {
    stop("Internal matching error: ncol(expr_matrix) != nrow(clinical_df)")
  }
  
  # Initialize results list
  results <- list()
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  # Standardize gender to "Male"/"Female" (ignore others/unknown)
  g <- tolower(trimws(as.character(clinical_df$gender)))
  gender_std <- rep(NA_character_, length(g))
  gender_std[grepl("^m", g)] <- "Male"
  gender_std[grepl("^f", g)] <- "Female"
  clinical_df$gender_std <- factor(gender_std, levels = c("Male", "Female"))
  
  # Loop through each gene
  for (gene in top15_genes_clean) {
    
    # Skip gene if not in the expression matrix
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    # Get gene symbol or fallback to Ensembl ID
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    # Create merged dataframe
    df <- clinical_df
    # attach expression vector (length == nrow(df) by construction)
    df$expression <- as.numeric(expr_matrix[gene, use_ids])
    
    # Check required survival columns exist
    if (!all(c("overall_survival", "deceased", "gender_std") %in% colnames(df))) {
      message("Missing required columns for gene: ", gene_name)
      next
    }
    
    # Filter out missing values
    df <- df[!is.na(df$expression) &
               !is.na(df$overall_survival) &
               !is.na(df$deceased) &
               !is.na(df$gender_std), , drop = FALSE]
    
    # Loop over gender groups
    for (sex in levels(df$gender_std)) {
      
      # Subset data for this gender
      sub_df <- df[df$gender_std == sex, , drop = FALSE]
      
      # Split into High/Low expression based on median
      sub_df$group <- factor(
        ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High")   # Low = reference in Cox
      )
      
      # Skip if too few samples or only one expression group
      if (nrow(sub_df) < 10 || length(unique(sub_df$group)) < 2) next
      
      # Fit Kaplan-Meier curve
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      
      # Log-rank test
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Cox regression
      cox <- tryCatch({
        summary(coxph(Surv(overall_survival, deceased) ~ group, data = sub_df))
      }, error = function(e) NULL)
      
      # If Cox regression ran successfully, store results
      if (!is.null(cox)) {
        results[[length(results) + 1]] <- data.frame(
          gene_id     = gene,
          gene_symbol = gene_name,
          gender      = sex,
          HR          = round(cox$coefficients[, "exp(coef)"], 3),
          lower_CI    = round(cox$conf.int[, "lower .95"], 3),
          upper_CI    = round(cox$conf.int[, "upper .95"], 3),
          KM_pval     = signif(pval, 3),
          Cox_pval    = signif(cox$coefficients[, "Pr(>|z|)"], 3),
          n           = nrow(sub_df)
        )
      }
      
      # Save plot if required
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        # filename-safe gender tag
        sex_tag <- tolower(gsub("[^A-Za-z0-9]+", "", sex))
        
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = paste(gene_name, "in", cancer_type, "- Gender", sex),
          xlab = "Days",
          ylab = "Survival Probability"
        )
        
        filename <- paste0(cancer_type, "_", gene_name, "_Gender_", sex_tag, "_survival.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        cat("Saved survival plot for:", gene_name, "- Gender", sex, "\n")
      }
    }
  }
  
  # Combine and export results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_gender_stratified_survival_results_ALL.csv"), row.names = FALSE)
    cat("Saved gender-stratified survival results for", nrow(results_df),
        "gene/gender combinations in", cancer_type, "\n")
  } else {
    cat("No valid survival results for", cancer_type, "\n")
  }
}

# Apply to LUAD
gender_stratified_survival_and_cox(expr_matrix = LUAD_expr,clinical_df = clinical_LUAD,top15_genes = top15_LUAD,
                                   master_gene_label = master_gene_label,cancer_type = "LUAD",save_all_plots = TRUE)

# Apply to LUSC
gender_stratified_survival_and_cox(expr_matrix = LUSC_expr,clinical_df = clinical_LUSC,top15_genes = top15_LUSC,
                                   master_gene_label = master_gene_label,cancer_type = "LUSC",save_all_plots = TRUE)


#------------------------- Smoking-Based Analysis of lncRNAs --------------------
# Load required libraries
library(dplyr)               # Data wrangling

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

# Perform smoking-status correlation for the 15 genes per cancer type
smoking_correlation_top15 <- function(expr_file, clinical_file, cancer_type, top15_genes, master_gene_label, sample_type = c("tumor", "normal")) {
  
  # Ensure sample_type is one of the allowed values
  sample_type <- match.arg(sample_type)
  
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1)
  clinical <- read.csv(clinical_file)
  
  # Filter expression columns based on sample type suffix (_tumor or _normal)
  if (sample_type == "tumor") {
    expr <- expr[, grepl("_tumor$", colnames(expr))]
  } else {
    expr <- expr[, grepl("_normal$", colnames(expr))]
  }
  
  # Fix sample IDs to match clinical$submitter_id
  fix_sample_ids <- function(ids) {
    sapply(ids, function(x) {
      parts <- unlist(strsplit(x, split = "\\."))
      parts <- parts[parts != ""]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
    })
  }
  
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Match clinical to expression samples
  clinical <- clinical[match(colnames(expr), clinical$submitter_id), ]
  
  # Recode smoking status into 3 groups: Never smoker, Ex-smoker, Current smoker
  clinical$smoking_group <- dplyr::case_when(
    clinical$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
    grepl("Reformed", clinical$tobacco_smoking_status) ~ "Ex-smoker",
    clinical$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
    TRUE ~ NA_character_
  )
  
  # Print counts per smoking group
  print(table(clinical$smoking_group))
  
  # Clean Top 15 gene IDs to remove version numbers
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Run Kruskal–Wallis test for each Top 15 gene
  results <- lapply(top15_genes_clean, function(gene) {
    if (gene %in% rownames(expr)) {
      expr_vals <- as.numeric(expr[gene, ])
      smoking_vals <- clinical$smoking_group
      
      # Only run test if at least 2 smoking groups are present
      if (length(unique(na.omit(smoking_vals))) >= 2) {
        test <- kruskal.test(expr_vals ~ smoking_vals)
        
        # Calculate absolute median difference and direction of change
        group_medians <- tapply(expr_vals, smoking_vals, median, na.rm = TRUE)
        abs_median_diff <- abs(max(group_medians, na.rm = TRUE) - min(group_medians, na.rm = TRUE))
        highest_group <- names(group_medians)[which.max(group_medians)]
        lowest_group <- names(group_medians)[which.min(group_medians)]
        direction <- paste0(highest_group, " > ", lowest_group)
        
        # Return results for this gene
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
  
  # Combine results into a data frame
  results_df <- do.call(rbind, results)
  
  # Check if results_df is NULL or empty
  if (is.null(results_df)) {
    cat("No valid smoking correlations computed for", cancer_type, "-", sample_type, "\n")
    return(NULL)
  }
  
  # Adjust p-values using Benjamini–Hochberg method
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Map Ensembl IDs to gene symbols if available
  results_df$gene_symbol <- ifelse(
    is.na(master_gene_label[results_df$gene_id]) | master_gene_label[results_df$gene_id] == "",
    results_df$gene_id,
    master_gene_label[results_df$gene_id]
  )
  
  # Save results
  filename <- paste0(cancer_type, "_", sample_type, "_smoking3_correlation_top15.csv")
  write.csv(results_df, filename, row.names = FALSE)
  cat("Done with", cancer_type, "-", sample_type, "- saved smoking correlation for Top 15 DEGs\n")
}

# Example usage:
# LUAD - tumor and normal
smoking_correlation_top15("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, sample_type = "tumor")
smoking_correlation_top15("LUAD_lncRNA_log_trans_after_norm.csv", "clinical_LUAD_clean.csv", "LUAD", top15_LUAD, master_gene_label, sample_type = "normal")

# LUSC - tumor and normal
smoking_correlation_top15("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, sample_type = "tumor")
smoking_correlation_top15("LUSC_lncRNA_log_trans_after_norm.csv", "clinical_LUSC_clean.csv", "LUSC", top15_LUSC, master_gene_label, sample_type = "normal")



#================= Plot Top 15 Smoking-Correlated DEGs =================
# This script generates boxplots for the Top 15 lncRNAs showing significant differences 
# in expression across smoking status groups ("Never smoker", "Ex-smoker", "Current smoker").
# The function works for either tumor or normal samples from a single expression file.

# Load required libraries
library(ggplot2)   # For plotting
library(ggpubr)    # For statistical annotations on plots
library(readr)     # For reading CSV files
library(dplyr)     # For data manipulation

# Function: fix_sample_ids
# Purpose: Ensure expression matrix sample names match TCGA clinical submitter_id format
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    # Split the string on '.' and remove empty entries
    parts <- unlist(strsplit(x, split = "\\."))
    parts <- parts[parts != ""]
    # Keep only the first 3 parts (TCGA-XX-XXXX) to match clinical data IDs
    if (length(parts) >= 3) paste(parts[1:3], collapse = "-") else NA
  })
}

# Function: add_smoking_group
# Purpose: Convert detailed smoking status into 3 consistent categories
add_smoking_group <- function(clinical_df) {
  clinical_df$smoking_group <- dplyr::case_when(
    clinical_df$tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
    grepl("Reformed", clinical_df$tobacco_smoking_status) ~ "Ex-smoker",
    clinical_df$tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
    TRUE ~ NA_character_
  )
  return(clinical_df)
}

# Function: plot_top15_smoking_genes
# Purpose: Generate and save boxplots for the Top 15 DEGs by smoking group
plot_top15_smoking_genes <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type, sample_type = c("tumor","normal")) {
  
  # Ensure sample_type is valid ("tumor" or "normal")
  sample_type <- match.arg(sample_type)
  
  # Load expression matrix
  expr <- read.csv(expr_file, row.names = 1)
  
  # Filter columns based on sample type suffix
  # e.g., "_tumor" for tumor samples, "_normal" for normal samples
  if (sample_type == "tumor") {
    expr <- expr[, grepl("_tumor$", colnames(expr)), drop = FALSE]
  } else {
    expr <- expr[, grepl("_normal$", colnames(expr)), drop = FALSE]
  }
  
  # Check if any columns remain after filtering
  if (ncol(expr) == 0) {
    message("No columns found for sample_type = ", sample_type, " in ", expr_file)
    return(invisible(NULL))
  }
  
  # Load clinical data
  clinical <- read.csv(clinical_file)
  
  # Standardize column names in expression matrix to match clinical IDs
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Add smoking group classification to clinical data
  clinical <- add_smoking_group(clinical)
  
  # Match samples between expression and clinical datasets
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  
  # Subset both datasets to matched samples
  expr <- expr[, common_samples, drop = FALSE]
  clinical <- clinical[match(common_samples, clinical$submitter_id), ]
  
  # Skip plotting if too few samples after matching
  if (length(common_samples) < 3) {
    message("Too few matched samples (", length(common_samples), ") for plotting: ", cancer_type, " - ", sample_type)
    return(invisible(NULL))
  }
  
  # Clean gene IDs (remove version numbers, e.g., ENSG00000.1 → ENSG00000)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Set consistent factor order for smoking groups
  clinical$smoking_group <- factor(
    clinical$smoking_group,
    levels = c("Never smoker", "Ex-smoker", "Current smoker")
  )
  
  # Generate plots for each Top 15 gene
  for (gene in top15_genes_clean) {
    if (gene %in% rownames(expr)) {
      
      # Create data frame with expression values and smoking group
      expr_vec <- as.numeric(expr[gene, ])
      df <- data.frame(expression = expr_vec, smoking_group = clinical$smoking_group)
      
      # Remove rows with missing smoking group info
      df <- df[!is.na(df$smoking_group), ]
      
      # Skip if insufficient data after filtering
      if (nrow(df) < 3 || length(unique(df$smoking_group)) < 2) {
        message("Skipping ", gene, ": insufficient groups after NA removal.")
        next
      }
      
      # Sample sizes for labeling above boxplots
      df %>%
        group_by(smoking_group) %>%
        summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.10, .groups = "drop") -> label_df
      
      # Only include pairwise comparisons for groups that are actually present in the data
      present_groups <- levels(droplevels(df$smoking_group))
      present_groups <- present_groups[!is.na(present_groups)]
      all_comparisons <- list(
        c("Never smoker", "Ex-smoker"),
        c("Never smoker", "Current smoker"),
        c("Ex-smoker", "Current smoker")
      )
      comparisons_list <- Filter(function(x) all(x %in% present_groups), all_comparisons)
      
      # Resolve gene symbol from master_gene_label mapping
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Build boxplot with statistical annotations
      p <- ggplot(df, aes(x = smoking_group, y = expression, fill = smoking_group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
        # Add Wilcoxon test p-values if we have valid group pairs
        { if (length(comparisons_list) > 0)
          stat_compare_means(
            comparisons = comparisons_list,
            method = "wilcox.test",
            label = "p.signif",
            bracket.size = 0.5,
            step.increase = 0.1,
            size = 5
          )
          else NULL } +
        # Add sample size labels above each box
        geom_text(
          data = label_df,
          aes(x = smoking_group, y = y, label = paste0("n = ", n)),
          inherit.aes = FALSE, size = 5.5, color = "black"
        ) +
        # Manual color scheme for consistency
        scale_fill_manual(values = c("Never smoker" = "grey", "Ex-smoker" = "#A80000", "Current smoker" = "#F8766D")) +
        coord_cartesian(ylim = c(NA, max(df$expression, na.rm = TRUE) * 1.4)) +
        labs(
          title = paste(gene_symbol, "expression in", cancer_type, "-", sample_type),
          x = "Smoking group",
          y = "Variance Stabilized Expression (VST)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 22),
          axis.title  = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title= element_text(size = 20),
          plot.title  = element_text(size = 22, face = "bold")
        )
      
      # Save plot to file
      ggsave(
        filename = paste0(cancer_type, "_", sample_type, "_", gene_symbol, "_Top15_Smoking_Boxplot.png"),
        plot = p, width = 10, height = 8, dpi = 300
      )
      
      cat("Saved plot for gene:", gene_symbol, "(", cancer_type, "-", sample_type, ")\n")
      
    } else {
      message("Gene not found in expression matrix: ", gene)
    }
  }
  
  cat("Finished plotting Top 15 smoking-correlated lncRNAs for", cancer_type, "-", sample_type, "\n")
}

# Example usage
plot_top15_smoking_genes(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,master_gene_label = master_gene_label,cancer_type = "LUAD",sample_type = "tumor")

plot_top15_smoking_genes(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD,master_gene_label = master_gene_label,cancer_type = "LUAD",sample_type = "normal")

plot_top15_smoking_genes(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUSC_clean.csv",top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,cancer_type = "LUSC",sample_type = "tumor")

plot_top15_smoking_genes(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv",clinical_file = "clinical_LUSC_clean.csv",top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,cancer_type = "LUSC",sample_type = "normal")



# ---------------- Smoking-Stratified Survival Analysis - Top 15 ---------------
# Smoking-Stratified Survival Analysis for Top 15 DE lncRNAs
# Load required libraries
library(survival)    # For survival models 
library(survminer)   # For visualization of survival curves with ggplot2
library(dplyr)       # For data manipulation
library(stringr)     # For safe filename text wrapping/cleaning (used in plot title)


#Load required data
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1, check.names = FALSE)
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1, check.names = FALSE)
clinical_LUAD <- read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)


# This function performs survival analysis within each smoking group
smoking_stratified_survival_and_cox <- function(expr_matrix, clinical_df, top15_genes,
                                                master_gene_label, cancer_type,
                                                save_all_plots = TRUE) {
  #ID HARMONIZATION & MATCHING (VERY IMPORTANT)
  # expr_matrix is genes x samples with TCGA-XX-XXXX as colnames.
  # clinical_df has submitter_id; we match on the first 12 chars and REORDER
  # the clinical rows to be in the SAME order as the expression columns.
  colnames(expr_matrix) <- gsub("\\.", "-", colnames(expr_matrix))                # Just in case
  clinical_df$submitter_id <- substr(as.character(clinical_df$submitter_id), 1, 12)
  
  use_ids <- intersect(colnames(expr_matrix), clinical_df$submitter_id)
  if (length(use_ids) < 10) {
    stop(paste("Too few matched samples for", cancer_type, ":", length(use_ids)))
  }
  
  expr_matrix <- expr_matrix[, use_ids, drop = FALSE]
  clinical_df <- clinical_df[match(use_ids, clinical_df$submitter_id), , drop = FALSE]
  
  if (ncol(expr_matrix) != nrow(clinical_df)) {
    stop("Internal matching error: ncol(expr_matrix) != nrow(clinical_df)")
  }
  
  # Initialize results list
  results <- list()
  
  # Clean gene IDs
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  # Ensure smoking status is available and create 3-group variable
  # (more robust recoding using case-insensitive matching)
  if (!"tobacco_smoking_status" %in% colnames(clinical_df)) {
    message("Missing 'tobacco_smoking_status' in clinical_df.")
    return(invisible(NULL))
  }
  smoking_raw <- tolower(trimws(as.character(clinical_df$tobacco_smoking_status)))
  clinical_df$smoking_group <- dplyr::case_when(
    grepl("non[- ]?smoker|never", smoking_raw)           ~ "Never smoker",
    grepl("reformed|former|ex", smoking_raw)             ~ "Ex-smoker",
    grepl("current", smoking_raw)                        ~ "Current smoker",
    TRUE                                                 ~ NA_character_
  )
  clinical_df$smoking_group <- factor(
    clinical_df$smoking_group,
    levels = c("Never smoker", "Ex-smoker", "Current smoker")
  )
  
  # Create/verify survival columns:
  # Expect either precomputed 'overall_survival' and 'deceased' OR columns to derive them.
  has_precomputed <- all(c("overall_survival", "deceased") %in% colnames(clinical_df))
  if (!has_precomputed) {
    needed <- c("vital_status", "days_to_death", "days_to_last_follow_up")
    if (!all(needed %in% colnames(clinical_df))) {
      message("Missing survival columns. Provide overall_survival & deceased, or vital_status/days_to_*.")
      return(invisible(NULL))
    }
    clinical_df$deceased <- as.numeric(clinical_df$vital_status != "Alive")  # 1 = event
    clinical_df$overall_survival <- ifelse(
      clinical_df$deceased == 1,
      as.numeric(clinical_df$days_to_death),
      as.numeric(clinical_df$days_to_last_follow_up)
    )
  }
  
  # Loop through each gene
  for (gene in top15_genes_clean) {
    
    # Skip gene if not in the expression matrix
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    # Get gene symbol or fallback to Ensembl ID
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    # Create merged dataframe
    df <- clinical_df
    # attach expression vector (length == nrow(df) by construction)
    df$expression <- as.numeric(expr_matrix[gene, use_ids])
    
    # Filter out missing values required for survival and smoking
    df <- df[!is.na(df$expression) &
               !is.na(df$overall_survival) &
               !is.na(df$deceased) &
               !is.na(df$smoking_group), , drop = FALSE]
    if (nrow(df) < 10) next
    
    # Loop over smoking groups
    for (smoking in levels(df$smoking_group)) {
      if (is.na(smoking)) next
      
      # Subset data for this smoking group
      sub_df <- df[df$smoking_group == smoking, , drop = FALSE]
      if (nrow(sub_df) < 10) next
      
      # Split into High/Low expression based on median
      sub_df$group <- factor(
        ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High")   # Low = reference in Cox
      )
      
      # Skip if only one group remains
      if (length(unique(sub_df$group)) < 2) next
      
      # Fit Kaplan-Meier curve
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      
      # Log-rank test
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Cox regression (High vs Low)
      cox <- tryCatch({
        summary(coxph(Surv(overall_survival, deceased) ~ group, data = sub_df))
      }, error = function(e) NULL)
      
      # Store results if Cox regression ran successfully
      if (!is.null(cox)) {
        results[[length(results) + 1]] <- data.frame(
          gene_id       = gene,
          gene_symbol   = gene_name,
          smoking_group = as.character(smoking),
          HR            = round(cox$coefficients[, "exp(coef)"], 3),
          lower_CI      = round(cox$conf.int[, "lower .95"], 3),
          upper_CI      = round(cox$conf.int[, "upper .95"], 3),
          KM_pval       = signif(pval, 3),
          Cox_pval      = signif(cox$coefficients[, "Pr(>|z|)"], 3),
          n             = nrow(sub_df)
        )
      } else {
        # Still record KM p-value if Cox fails
        results[[length(results) + 1]] <- data.frame(
          gene_id       = gene,
          gene_symbol   = gene_name,
          smoking_group = as.character(smoking),
          HR            = NA,
          lower_CI      = NA,
          upper_CI      = NA,
          KM_pval       = signif(pval, 3),
          Cox_pval      = NA,
          n             = nrow(sub_df)
        )
      }
      
      # Save plot if required (all or only significant)
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        # optional: title wrapping (use stringr if available)
        title_txt <- paste0(gene_name, " in ", cancer_type, " - Smoking: ", smoking)
        if (requireNamespace("stringr", quietly = TRUE)) {
          title_txt <- stringr::str_wrap(title_txt, width = 50)
        }
        
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = title_txt,
          xlab = "Days",
          ylab = "Survival Probability"
        )
        safe_smoking <- gsub("[^A-Za-z0-9]+", "_", as.character(smoking))
        filename <- paste0(cancer_type, "_", gene_name, "_Smoking_", safe_smoking, "_survival.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        cat("Saved survival plot for:", gene_name, "- Smoking:", smoking, "\n")
      }
    }
  }
  
  # Combine and export results
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_smoking_stratified_survival_results_ALL.csv"), row.names = FALSE)
    cat("Saved smoking-stratified survival results for", nrow(results_df),
        "gene/smoking group combinations in", cancer_type, "\n")
  } else {
    cat("No valid smoking-stratified survival results for", cancer_type, "\n")
  }
}

# Apply to LUAD
smoking_stratified_survival_and_cox(
  expr_matrix = LUAD_expr,clinical_df = clinical_LUAD,top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,cancer_type = "LUAD",save_all_plots = TRUE)

# Apply to LUSC
smoking_stratified_survival_and_cox(
  expr_matrix = LUSC_expr,clinical_df = clinical_LUSC,top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,cancer_type = "LUSC",save_all_plots = TRUE)

# ------------------------ Stage-Based Analysis of lncRNAs ------------------------
# Kruskal-Wallis test to identify DE lncRNAs significantly different by stage (I, II, III, IV)
# Stage-Based Analysis of lncRNAs (Top 15) with Kruskal–Wallis tests
# Goal: Test whether expression of the Top 15 lncRNAs differs by AJCC pathologic stage (I–IV)

# Load necessary libraries
library(dplyr)   # Data wrangling
library(readr)   # Reading CSV files

# Function: fix_sample_ids
# Purpose: Convert expression column names to TCGA submitter_id format (TCGA-XX-XXXX)
fix_sample_ids <- function(ids) {
  ids <- gsub("(_tumor|_normal)$", "", ids, ignore.case = TRUE)  # drop suffix if present
  ids <- gsub("\\.", "-", ids)                                   # dots -> dashes
  vapply(strsplit(ids, "-"), function(p) {
    p <- p[p != ""]
    if (length(p) >= 3) paste(p[1:3], collapse = "-") else NA_character_
  }, character(1))
}

# Helper: normalize AJCC pathologic stage to "Stage I/II/III/IV" (drops sub-stages IA/IB, etc.)
normalize_stage <- function(x) {
  x <- as.character(x)
  # pull the roman numeral after 'Stage' (case-insensitive), ignore trailing letters/numbers
  roman <- toupper(sub(".*[Ss]tage\\s*([IVX]+).*", "\\1", x))
  # keep only I/II/III/IV; everything else -> NA
  roman[!roman %in% c("I","II","III","IV")] <- NA
  out <- ifelse(is.na(roman), NA, paste("Stage", roman))
  factor(out, levels = c("Stage I","Stage II","Stage III","Stage IV"))
}

# Function: stage_correlation_top15
# Purpose: For each of the Top 15 lncRNAs, run a Kruskal–Wallis test across pathologic stages.
# Inputs: expr_file, clinical_file, top15_genes, master_gene_label, cancer_type
stage_correlation_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression data and clinical data
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  
  # Standardize expression column names to TCGA format (TCGA-XX-XXXX)
  colnames(expr) <- fix_sample_ids(colnames(expr))
  # Standardize clinical submitter_id to first 12 chars and same punctuation
  clinical$submitter_id <- substr(gsub("\\.", "-", clinical$submitter_id), 1, 12)
  
  # Align expression and clinical samples (intersection ensures matched order later)
  shared_samples <- intersect(colnames(expr), clinical$submitter_id)
  if (length(shared_samples) < 10) {
    stop(paste0("[", cancer_type, "] Too few matched samples: ", length(shared_samples)))
  }
  expr <- expr[, shared_samples, drop = FALSE]
  clinical <- clinical[match(shared_samples, clinical$submitter_id), , drop = FALSE]
  
  # Normalize/clean AJCC stage labels into Stage I/II/III/IV (drop sub-stages like IA/IB)
  clinical$stage_group <- normalize_stage(clinical$ajcc_pathologic_stage)
  
  # Report sample counts per stage (useful sanity check)
  message("[", cancer_type, "] samples per stage (including NA):")
  print(table(clinical$stage_group, useNA = "ifany"))
  
  # Clean gene IDs: remove version suffix from Top 15 and rownames
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # For each gene, run Kruskal–Wallis across available stages (requires >= 2 groups present)
  results <- lapply(top15_genes_clean, function(gene) {
    if (!gene %in% rownames(expr)) {
      message("Gene not found in expression matrix: ", gene)
      return(NULL)
    }
    
    # pairwise-complete: drop rows with NA stage or NA expression
    expr_vals  <- as.numeric(expr[gene, ])
    stage_vals <- clinical$stage_group
    ok <- !is.na(expr_vals) & !is.na(stage_vals)
    expr_vals  <- expr_vals[ok]
    stage_vals <- droplevels(stage_vals[ok])
    
    # Only proceed if at least two stages are represented
    if (length(unique(stage_vals)) < 2) return(NULL)
    
    # Kruskal–Wallis test
    test <- kruskal.test(expr_vals ~ stage_vals)
    
    # Compute absolute difference between max and min stage medians (effect size proxy)
    group_medians <- tapply(expr_vals, stage_vals, median, na.rm = TRUE)
    abs_median_diff <- abs(max(group_medians, na.rm = TRUE) - min(group_medians, na.rm = TRUE))
    
    # Map to gene symbol if available, else keep Ensembl ID
    gene_symbol <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    data.frame(
      gene_id         = gene,
      gene_symbol     = gene_symbol,
      p_value         = unname(test$p.value),
      abs_median_diff = unname(abs_median_diff),
      n_total         = length(expr_vals)
    )
  })
  
  # Combine results across genes
  results_df <- do.call(rbind, results)
  
  # Exit cleanly if nothing to report
  if (is.null(results_df) || nrow(results_df) == 0) {
    cat("No valid results to process (insufficient stage variation or unmatched samples)\n")
    return(NULL)
  }
  
  # Multiple testing correction (Benjamini–Hochberg)
  results_df$padj <- p.adjust(results_df$p_value, method = "BH")
  
  # Save all results and significant subset
  write.csv(results_df, paste0(cancer_type, "_stage_correlation_top15_all.csv"), row.names = FALSE)
  write.csv(results_df[results_df$padj < 0.05, ], paste0(cancer_type, "_stage_correlation_top15_significant.csv"), row.names = FALSE)
  
  cat("Done with", cancer_type, "- Top 15 stage correlation:",
      sum(results_df$padj < 0.05), "significant genes\n")
}

# Example calls
stage_correlation_top15("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD")
stage_correlation_top15("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC")

#=============Plot stage-wise correlation
# Plot Top 15 lncRNAs by Pathologic Stage (boxplots + pairwise Wilcoxon)
# Goal: Visualize expression differences across Stage I/II/III/IV for each Top 15 lncRNA

# Required libraries
library(ggplot2)  # Plotting
library(ggpubr)   # stat_compare_means for pairwise tests
library(readr)    # Reading CSV files
library(dplyr)    # Data manipulation

# Function: fix_sample_ids
# Purpose: Convert expression column names to TCGA submitter_id format (TCGA-XX-XXXX)
fix_sample_ids <- function(ids) {
  ids <- gsub("(_tumor|_normal)$", "", ids, ignore.case = TRUE)  # drop suffix if present
  ids <- gsub("\\.", "-", ids)                                   # dots -> dashes
  vapply(strsplit(ids, "-"), function(p) {
    p <- p[p != ""]
    if (length(p) >= 3) paste(p[1:3], collapse = "-") else NA_character_
  }, character(1))
}

# Helper: collapse duplicate columns (same submitter_id) by mean
collapse_dupes_by_mean <- function(mat) {
  ids <- colnames(mat)
  idx <- split(seq_along(ids), ids)
  out <- lapply(idx, function(j)
    if (length(j) == 1) as.numeric(mat[, j])
    else rowMeans(as.matrix(mat[, j, drop = FALSE]), na.rm = TRUE)
  )
  out <- do.call(cbind, out)
  colnames(out) <- names(idx); rownames(out) <- rownames(mat)
  out
}

# normalize AJCC pathologic stage to "Stage I/II/III/IV" (drops IA/IB, etc.)
normalize_stage <- function(x) {
  x <- as.character(x)
  roman <- toupper(sub(".*[Ss]tage\\s*([IVX]+).*", "\\1", x))
  roman[!roman %in% c("I","II","III","IV")] <- NA
  out <- ifelse(is.na(roman), NA, paste("Stage", roman))
  factor(out, levels = c("Stage I","Stage II","Stage III","Stage IV"))
}

# Function: plot_top_stage_genes_top15
# Purpose: Produce a boxplot per gene showing expression across Stage I–IV with pairwise Wilcoxon p-values.
plot_top_stage_genes_top15 <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  
  # Load expression data and clinical metadata
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  
  #keep tumor-only columns if the file contains both tumor/normal ---
  tumor_cols <- grep("_tumor$", colnames(expr), value = TRUE, ignore.case = TRUE)
  if (length(tumor_cols) > 0) {
    expr <- expr[, tumor_cols, drop = FALSE]
  }
  
  # Standardize expression column names to TCGA format and match to clinical
  colnames(expr) <- fix_sample_ids(colnames(expr))
  expr <- expr[, !is.na(colnames(expr)), drop = FALSE]
  # Collapse duplicate submitter IDs (replicates) by mean
  expr <- collapse_dupes_by_mean(expr)
  
  # Clinical IDs: ensure same 12-char submitter_id format
  clinical$submitter_id <- substr(gsub("\\.", "-", clinical$submitter_id), 1, 12)
  
  # Normalize/clean AJCC stage labels into Stage I/II/III/IV
  clinical$stage_group <- normalize_stage(clinical$ajcc_pathologic_stage)
  
  # Align on shared samples (keeps order consistent in both)
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  if (length(common_samples) < 10) {
    stop(paste0("[", cancer_type, "] Too few matched tumor samples: ", length(common_samples)))
  }
  expr <- expr[, common_samples, drop = FALSE]
  clinical <- clinical[match(common_samples, clinical$submitter_id), , drop = FALSE]
  
  # Clean gene IDs (remove version numbers from Top 15 and rownames)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Loop through Top 15 genes and build plots
  for (gene in top15_genes_clean) {
    
    if (gene %in% rownames(expr)) {
      # Build plotting dataframe: expression + stage_group
      expr_vec <- as.numeric(expr[gene, ])
      df <- data.frame(
        expression = expr_vec,
        stage_group = clinical$stage_group
      )
      # Drop samples without stage info or expression
      df <- df[!is.na(df$stage_group) & !is.na(df$expression), , drop = FALSE]
      # Skip plotting if fewer than 2 distinct stages remain
      if (nrow(df) < 6 || length(unique(df$stage_group)) < 2) next
      
      # Sample sizes and label positions per stage (for n annotations)
      label_df <- df %>%
        group_by(stage_group) %>%
        summarise(n = n(), y = max(expression, na.rm = TRUE) * 1.12, .groups = "drop")
      
      # Define pairwise comparisons among stages; keep only present pairs
      all_comparisons <- list(
        c("Stage I", "Stage II"),
        c("Stage I", "Stage III"),
        c("Stage I", "Stage IV"),
        c("Stage II", "Stage III"),
        c("Stage II", "Stage IV"),
        c("Stage III", "Stage IV")
      )
      present_levels <- levels(droplevels(df$stage_group))
      comparisons_list <- Filter(function(x) all(x %in% present_levels), all_comparisons)
      
      # Resolve gene symbol for plotting
      gene_symbol <- ifelse(
        is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
        gene,
        master_gene_label[gene]
      )
      
      # Build boxplot with Wilcoxon pairwise p-values and sample sizes
      p <- ggplot(df, aes(x = stage_group, y = expression, fill = stage_group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        { if (length(comparisons_list) > 0)
          stat_compare_means(
            comparisons = comparisons_list,
            method = "wilcox.test",
            label = "p.signif",
            tip.length = 0.01,
            step.increase = 0.1,
            size = 5
          )
          else NULL } +
        geom_text(
          data = label_df,
          aes(x = stage_group, y = y, label = paste0("n = ", n)),
          inherit.aes = FALSE, size = 5.5, color = "black"
        ) +
        scale_fill_manual(values = c(
          "Stage I"   = "red",
          "Stage II"  = "grey",
          "Stage III" = "#A80000",
          "Stage IV"  = "#D62728"
        )) +
        labs(
          title = paste(gene_symbol, "expression by stage in", cancer_type),
          x = "Pathologic Stage",
          y = "Variance Stabilized Expression (VST)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text  = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 14),
          plot.title   = element_text(size = 18, face = "bold")
        )
      
      # Save plot to disk (one PNG per gene)
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

# Example calls
plot_top_stage_genes_top15(
  expr_file = "LUAD_lncRNA_log_trans_after_norm.csv", clinical_file = "clinical_LUAD_clean.csv",
  top15_genes = top15_LUAD, master_gene_label = master_gene_label, cancer_type = "LUAD"
)

plot_top_stage_genes_top15(
  expr_file = "LUSC_lncRNA_log_trans_after_norm.csv", clinical_file = "clinical_LUSC_clean.csv",
  top15_genes = top15_LUSC, master_gene_label = master_gene_label, cancer_type = "LUSC"
)


# ------------------------ Stage-Based Survival Analysis of lncRNAs ------------------------
# Survival Analysis for DE lncRNAs by Pathologic Stage
# Load required libraries
library(survival)    # For survival models (Kaplan-Meier, Cox regression)
library(survminer)   # For visualization of survival curves with ggplot2
library(dplyr)       # For data manipulation
library(stringr)     # For safe plot titles / filenames

#Load datasets
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1, check.names = FALSE)
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1, check.names = FALSE)

clinical_LUAD <- readr::read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- readr::read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)

# This function:
# Takes the top 15 lncRNAs for a given cancer type
# For each gene, splits patients by AJCC pathologic stage (I–IV)
# Within each stage, separates samples into "High" and "Low" expression groups (median split)
# Runs Kaplan–Meier survival analysis and log-rank test and Cox proportional hazards regression to calculate hazard ratios (HR) and confidence intervals (CI)

stage_stratified_survival_and_cox <- function(expr_matrix, clinical_df, top15_genes,
                                              master_gene_label, cancer_type,
                                              save_all_plots = TRUE) {
  
  # Initialize empty list to store results
  results <- list()
  
  # Remove version numbers from Ensembl IDs in Top 15 list (e.g., ENSG000001.1 -> ENSG000001)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  # Remove version numbers from expression matrix rownames
  rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))
  
  # IMPORTANT: harmonize IDs & match samples before any analysis
  # clinical submitter_id -> first 12 chars (TCGA-XX-XXXX); also normalize dots to dashes
  clin_ids <- substr(gsub("\\.", "-", clinical_df$submitter_id), 1, 12)
  # we assume expr_matrix columns are already TCGA-XX-XXXX (from your *_expr_matched.csv files)
  use_ids <- intersect(colnames(expr_matrix), clin_ids)
  if (length(use_ids) < 10) {
    stop(paste0("[", cancer_type, "] Too few matched samples after alignment: ", length(use_ids)))
  }
  # align order in both objects
  expr_matrix <- expr_matrix[, use_ids, drop = FALSE]
  clinical_df <- clinical_df[match(use_ids, clin_ids), , drop = FALSE]
  
  # Sanity check
  if (ncol(expr_matrix) != nrow(clinical_df)) {
    stop("Internal matching error: ncol(expr_matrix) != nrow(clinical_df)")
  }
  
  # Ensure survival columns exist
  # If 'overall_survival' and 'deceased' columns are missing, try creating them from available clinical data
  has_required <- all(c("overall_survival", "deceased") %in% colnames(clinical_df))
  if (!has_required) {
    needed <- c("vital_status", "days_to_death", "days_to_last_follow_up")
    if (!all(needed %in% colnames(clinical_df))) {
      message("Missing survival columns: provide 'overall_survival' & 'deceased', or 'vital_status' + days_to_* fields.")
      return(invisible(NULL))
    }
    # deceased: 1 if vital_status != "Alive", else 0
    clinical_df$deceased <- as.numeric(clinical_df$vital_status != "Alive")
    # overall_survival: use days_to_death if deceased, otherwise days_to_last_follow_up
    clinical_df$overall_survival <- ifelse(
      clinical_df$deceased == 1,
      as.numeric(clinical_df$days_to_death),
      as.numeric(clinical_df$days_to_last_follow_up)
    )
  }
  # Make sure survival columns are numeric
  clinical_df$overall_survival <- as.numeric(clinical_df$overall_survival)
  clinical_df$deceased         <- as.numeric(clinical_df$deceased)
  
  # Standardize stage labels
  # Convert stage values to a clean, consistent format: Stage I, Stage II, Stage III, Stage IV
  normalize_stage <- function(x) {
    x <- as.character(x)
    roman <- toupper(sub(".*[Ss]tage\\s*([IVX]+).*", "\\1", x))
    roman[!roman %in% c("I","II","III","IV")] <- NA
    out <- ifelse(is.na(roman), NA, paste("Stage", roman))
    factor(out, levels = c("Stage I","Stage II","Stage III","Stage IV"))
  }
  clinical_df$stage_group <- normalize_stage(clinical_df$ajcc_pathologic_stage)
  
  #For each gene in the Top 15, run survival analysis within each stage
  for (gene in top15_genes_clean) {
    
    # Skip if gene is missing from expression matrix
    if (!(gene %in% rownames(expr_matrix))) {
      message("Gene not found in expression matrix: ", gene)
      next
    }
    
    # Retrieve gene symbol from annotation table or fallback to Ensembl ID
    gene_name <- ifelse(
      is.na(master_gene_label[gene]) | master_gene_label[gene] == "",
      gene,
      master_gene_label[gene]
    )
    
    # Merge clinical and expression data (ordering already aligned)
    df <- clinical_df
    df$expression <- as.numeric(expr_matrix[gene, use_ids])
    
    # Remove rows with missing data in expression, survival, or stage
    df <- df[!is.na(df$expression) & !is.na(df$overall_survival) & !is.na(df$deceased) & !is.na(df$stage_group), ,
             drop = FALSE]
    
    # Skip if too few samples overall
    if (nrow(df) < 10) next
    
    #WITHIN-STAGE ANALYSES
    # Loop over stages actually present
    for (stg in levels(droplevels(df$stage_group))) {
      sub_df <- df[df$stage_group == stg, , drop = FALSE]
      if (nrow(sub_df) < 10) next  # Skip if too few samples in this stage
      
      # Median split into High and Low expression groups (within stage)
      sub_df$group <- ifelse(sub_df$expression >= median(sub_df$expression, na.rm = TRUE), "High", "Low")
      sub_df$group <- factor(sub_df$group, levels = c("Low","High"))  # set reference to "Low"
      
      # Skip if only one expression group present
      if (length(unique(sub_df$group)) < 2) next
      
      # Kaplan–Meier survival analysis
      fit <- survfit(Surv(overall_survival, deceased) ~ group, data = sub_df)
      
      # Log-rank test for significance
      test <- survdiff(Surv(overall_survival, deceased) ~ group, data = sub_df)
      pval <- 1 - pchisq(test$chisq, df = 1)
      
      # Cox proportional hazards regression
      cox <- tryCatch(
        summary(coxph(Surv(overall_survival, deceased) ~ group, data = sub_df)),
        error = function(e) NULL
      )
      
      # Store results (prefer Cox if available; still record KM p if not)
      if (!is.null(cox)) {
        results[[length(results) + 1]] <- data.frame(
          gene_id      = gene,
          gene_symbol  = gene_name,
          stage_group  = as.character(stg),
          HR           = round(cox$coefficients[, "exp(coef)"], 3),
          lower_CI     = round(cox$conf.int[, "lower .95"], 3),
          upper_CI     = round(cox$conf.int[, "upper .95"], 3),
          KM_pval      = signif(pval, 3),
          Cox_pval     = signif(cox$coefficients[, "Pr(>|z|)"], 3),
          n            = nrow(sub_df)
        )
      } else {
        results[[length(results) + 1]] <- data.frame(
          gene_id      = gene,
          gene_symbol  = gene_name,
          stage_group  = as.character(stg),
          HR           = NA,
          lower_CI     = NA,
          upper_CI     = NA,
          KM_pval      = signif(pval, 3),
          Cox_pval     = NA,
          n            = nrow(sub_df)
        )
      }
      
      # Save plot if required
      if (save_all_plots || (!is.na(pval) && pval < 0.05)) {
        plot <- ggsurvplot(
          fit,
          data = sub_df,
          pval = TRUE,
          title = paste(gene_name, "in", cancer_type, "-", stg),
          xlab = "Days",
          ylab = "Survival Probability"
        )
        # Safe filename
        safe_stage <- gsub("[^A-Za-z0-9]+", "_", as.character(stg))
        filename <- paste0(cancer_type, "_", gene_name, "_Stage_", safe_stage, "_survival.png")
        ggsave(filename, plot = plot$plot, width = 6, height = 5, dpi = 300)
        cat("Saved survival plot for:", gene_name, "-", stg, "\n")
      }
    }
  }
  
  #EXPORT
  if (length(results) > 0) {
    results_df <- do.call(rbind, results)
    write.csv(results_df, paste0(cancer_type, "_stage_stratified_survival_results_ALL.csv"), row.names = FALSE)
    cat("Saved stage-stratified survival results for", nrow(results_df), "gene/stage combinations in", cancer_type, "\n")
  } else {
    cat("No valid stage-stratified survival results for", cancer_type, "\n")
  }
}


# Apply function to LUAD
stage_stratified_survival_and_cox(
  expr_matrix = LUAD_expr,clinical_df = clinical_LUAD,top15_genes = top15_LUAD,
  master_gene_label = master_gene_label,cancer_type = "LUAD",save_all_plots = TRUE)

# Apply function to LUSC
stage_stratified_survival_and_cox(
  expr_matrix = LUSC_expr,clinical_df = clinical_LUSC,top15_genes = top15_LUSC,
  master_gene_label = master_gene_label,cancer_type = "LUSC",save_all_plots = TRUE)



#---------------------------- Multivariate Cox Regression ----------------------
# Multivariate Cox Regression for Top 15 lncRNAs (adjusting for age, gender, smoking, stage)
# Load required libraries
library(survival)  # Cox models and survival objects
library(dplyr)     # Data wrangling
library(readr)     # CSV I/O
library(stringr)   # String utilities (pattern matching)
library(forcats)   # Factor handling (if needed later)


#Load datasets
LUAD_expr <- read.csv("LUAD_expr_matched.csv", row.names = 1, check.names = FALSE)
LUSC_expr <- read.csv("LUSC_expr_matched.csv", row.names = 1, check.names = FALSE)

clinical_LUAD <- readr::read_csv("clinical_LUAD_clean.csv", show_col_types = FALSE)
clinical_LUSC <- readr::read_csv("clinical_LUSC_clean.csv", show_col_types = FALSE)

#Helper: standardize expression column sample IDs to TCGA submitter_id format
# - Replace dots with dashes
# - Keep the first 12 characters (TCGA-XX-XXXX)
fix_sample_ids <- function(ids) {
  ids <- gsub("\\.", "-", ids)   # Replace '.' with '-'
  ids <- substr(ids, 1, 12)      # Truncate to TCGA-XX-XXXX
  return(ids)
}

# Multivariate Cox regression function
#Inputs: expr_file, clinical_file, top15_genes, master_gene_label, cancer_type
run_multivariate_cox <- function(expr_file, clinical_file, top15_genes, master_gene_label, cancer_type) {
  # Load expression and clinical data
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  
  # Clean IDs
  #Expression colnames match TCGA submitter_id format
  # (robust: handle both already-matched IDs and raw IDs with dots)
  expr_ids <- gsub("\\.", "-", colnames(expr))
  expr_ids <- vapply(strsplit(expr_ids, "-"), function(p) {
    if (length(p) >= 3) paste(p[1:3], collapse = "-") else NA_character_
  }, character(1))
  colnames(expr) <- expr_ids
  
  #Drop version numbers from gene IDs in expression matrix
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  # Clean Top 15 gene IDs (remove version numbers)
  top15_genes_clean <- gsub("\\..*", "", top15_genes)
  
  # Align samples present in both expression and clinical data and keeps ordering consistent
  clinical$submitter_id <- substr(gsub("\\.", "-", as.character(clinical$submitter_id)), 1, 12)
  common_samples <- intersect(colnames(expr), clinical$submitter_id)
  expr <- expr[, common_samples, drop = FALSE]
  clinical <- clinical[match(common_samples, clinical$submitter_id), , drop = FALSE]
  
  # Bail out early if too few samples
  if (length(common_samples) < 20) {
    message("[", cancer_type, "] Too few matched samples for multivariate Cox: ", length(common_samples))
    return(invisible(NULL))
  }
  
  # Build/clean clinical covariates required for the multivariable model
  # - Convert age to years
  # - Set factors for gender, smoking group, and stage
  # - Create survival time and event indicator
  suppressWarnings({
    clinical <- dplyr::mutate(
      clinical,
      age = as.numeric(age_at_diagnosis) / 365,  # Age in years
      gender = factor(gender),                   # Ensure categorical
      # Collapse detailed smoking into 3 groups
      smoking_group = dplyr::case_when(
        tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never smoker",
        grepl("Reformed", tobacco_smoking_status %||% "", ignore.case = TRUE) ~ "Ex-smoker",
        tobacco_smoking_status == "Current Smoker" ~ "Current smoker",
        TRUE ~ NA_character_
      ),
      smoking_group = factor(smoking_group, levels = c("Never smoker","Ex-smoker","Current smoker")),
      # Normalize AJCC pathologic stage to Stage I/II/III/IV
      stage_group = gsub("^Stage\\s*([IV]+)[A-Z]?$", "Stage \\1", ajcc_pathologic_stage %||% ""),
      stage_group = factor(stage_group, levels = c("Stage I","Stage II","Stage III","Stage IV")),
      # Survival outcome: event and time (in days)
      deceased = as.numeric(vital_status != "Alive"),
      overall_survival = ifelse(
        deceased == 1, suppressWarnings(as.numeric(days_to_death)),
        suppressWarnings(as.numeric(days_to_last_follow_up))
      )
    )
  })
  
  # Container for per-gene Cox summaries
  results_list <- list()
  
  # Iterate over each target gene
  for (gene in top15_genes_clean) {
    # Skip if gene not present in expression matrix
    if (!(gene %in% rownames(expr))) {
      cat("Skipping", gene, "- not found in expression matrix\n")
      next
    }
    
    # Add gene expression as a column to clinical for modeling
    clinical$expression <- as.numeric(expr[gene, ])
    
    # Keep only complete cases across all required variables
    df <- dplyr::filter(
      clinical,
      !is.na(expression),
      !is.na(age),
      !is.na(gender),
      !is.na(smoking_group),
      !is.na(stage_group),
      !is.na(overall_survival),
      !is.na(deceased)
    )
    df <- droplevels(df)
    
    cat(" Gene:", gene, "- Valid samples:", nrow(df), "\n")
    
    # Require a minimum sample size for stability (rule-of-thumb threshold)
    if (nrow(df) < 20) {
      cat(" Skipping", gene, "- not enough valid samples (<20)\n")
      next
    }
    
    # Fit multivariable Cox model:
    # Surv(time, event) ~ expression + age + gender + smoking_group + stage_group
    # - expression is the effect of interest adjusted for covariates
    cox_model <- tryCatch({
      survival::coxph(survival::Surv(overall_survival, deceased) ~
                        expression + age + gender + smoking_group + stage_group,
                      data = df)
    }, error = function(e) {
      message("Failed model for gene: ", gene, " - ", e$message)
      return(NULL)
    })
    
    # If model fit succeeded, extract coefficient table
    if (!is.null(cox_model)) {
      # Map Ensembl -> symbol; fall back to Ensembl if missing
      gene_symbol <- master_gene_label[gene]
      if (is.null(gene_symbol) || is.na(gene_symbol) || gene_symbol == "") gene_symbol <- gene
      
      summary_cox <- summary(cox_model)
      coefs <- as.data.frame(summary_cox$coefficients)  # includes coef, exp(coef), se(coef), z, Pr(>|z|)
      coefs$gene_id <- gene
      coefs$gene_symbol <- gene_symbol
      coefs$variable <- rownames(coefs)  # which term (expression, age, gender..., stage...)
      rownames(coefs) <- NULL
      results_list[[gene]] <- coefs
    }
  }
  
  # If any gene succeeded, bind and write results
  if (length(results_list) > 0) {
    all_results <- dplyr::bind_rows(results_list)
    
    # Focus summary on the expression effect for per-gene significance reporting
    expression_rows <- all_results %>%
      dplyr::filter(variable == "expression") %>%
      dplyr::mutate(
        HR = exp(coef),                                           # hazard ratio for expression
        lower_CI = exp(coef - 1.96 * `se(coef)`),                 # Wald 95% CI lower
        upper_CI = exp(coef + 1.96 * `se(coef)`),                 # Wald 95% CI upper
        is_significant = `Pr(>|z|)` < 0.05 & (lower_CI > 1 | upper_CI < 1) # p<0.05 and CI not crossing 1
      )
    
    # Count significant genes (expression term only)
    n_sig <- sum(expression_rows$is_significant, na.rm = TRUE)
    
    # Save full coefficient tables (all variables) for transparency
    write.csv(all_results, paste0(cancer_type, "_multivariate_cox_results_top15.csv"), row.names = FALSE)
    
    # Console summary
    if (n_sig > 0) {
      cat(" Found", n_sig, "significant lncRNAs for", cancer_type, "\n")
    } else {
      cat(" No significant lncRNAs found for", cancer_type, "\n")
    }
  }
}


# Run for LUAD and LUSC
run_multivariate_cox("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv", top15_LUAD, master_gene_label, "LUAD")
run_multivariate_cox("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv", top15_LUSC, master_gene_label, "LUSC")


#---------------------------- Forest Plot ---------------------------------------
# Purpose: Visualize adjusted hazard ratios (HR) from the multivariate Cox models
#          for the Top 15 lncRNAs in LUAD and LUSC, with 95% CIs and significance.

# Load required libraries
library(readr)    # read_csv
library(dplyr)    # bind_rows, filter, mutate, arrange
library(forcats)  # fct_reorder for ordered factors
library(ggplot2)  # plotting

# More precise printing (helpful for debugging very close CI bounds)
options(digits = 10)

# Load multivariate Cox results for LUAD and LUSC
df_luad <- read_csv("LUAD_multivariate_cox_results_top15.csv")
df_lusc <- read_csv("LUSC_multivariate_cox_results_top15.csv")

# Tag each dataset with its cancer type
df_luad$cancer <- "LUAD"
df_lusc$cancer <- "LUSC"

# Combine and keep only the 'expression' rows (gene effect of interest)
# Convert log-HR (coef) + SE to HR and 95% CI on the HR scale
#   Build a display label and a significance flag
combined_df <- bind_rows(df_luad, df_lusc) %>%
  filter(variable == "expression") %>%                       # focus on gene expression term
  mutate(
    HR = exp(coef),                                          # transform log-HR to HR
    lower_CI = round(exp(coef - 1.96 * `se(coef)`), 3),      # Wald 95% CI lower bound
    upper_CI = round(exp(coef + 1.96 * `se(coef)`), 3),      # Wald 95% CI upper bound
    gene_display = fct_reorder(                              # order genes by HR (ascending)
      paste0(gene_symbol, " (", cancer, ")"), HR
    ),
    significance = ifelse(                                   # flag significant if p<0.05 and CI excludes 1
      `Pr(>|z|)` < 0.05 & (lower_CI > 1 | upper_CI < 1),
      "Significant", "Not significant"
    )
  )

# Quick debug/inspection printout of key values (optional, helpful for QC)
print(
  combined_df %>%
    select(gene_symbol, cancer, HR, lower_CI, upper_CI, `Pr(>|z|)`, significance) %>%
    arrange(significance, desc(HR)),
  n = 30
)

# Compute dynamic x-axis limits so all CIs fit nicely with a small margin
x_min <- min(combined_df$lower_CI) - 0.005
x_max <- max(combined_df$upper_CI) + 0.005

# Build the forest plot:
#    - point = HR
#    - horizontal error bar = 95% CI
#    - vertical dashed line at HR=1 (no effect)
#    - color by significance
#    - coordinate flip for classic forest-plot look (labels on y-axis)
p <- ggplot(combined_df, aes(x = HR, y = gene_display, color = significance)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_continuous(                                       # note: this is linear scale (see tip below)
    limits = c(x_min, x_max),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_color_manual(values = c("Not significant" = "gray60", "Significant" = "#A80000")) +
  labs(
    title = "Multivariate Cox Model (Top 15 lncRNAs)",
    x = "Hazard Ratio",                                     # label says HR; remove "(log scale)" since not log
    y = NULL,
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

# Save the figure
ggsave("multivariate_forest_plot_FIXED.png", plot = p, width = 10, height = 10, dpi = 300)



#----------------- Evaluate 2- and 3-gene lncRNA signatures for survival -------

# Evaluate 2- and 3-gene lncRNA signatures for survival (using precomputed OS/event)

# Required libraries
library(survival)
library(survminer)
library(utils)
library(dplyr)
library(readr)
library(stringr)

# Fix TCGA sample IDs: replace dots with hyphens & truncate to TCGA-XX-XXXX
fix_sample_ids <- function(ids) {
  sapply(ids, function(x) {
    id <- gsub("\\.", "-", x)
    substr(id, 1, 12)
  }, USE.NAMES = FALSE)
}

# Build a survival-ready dataframe using overall_survival and deceased
prepare_survival_df_precomputed <- function(expr_file, clin_file) {
  # Load data
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  clin <- read_csv(clin_file, show_col_types = FALSE)
  
  # Strip version numbers from Ensembl IDs so column names (after transpose) are clean
  rownames(expr) <- gsub("\\..*", "", rownames(expr))
  
  # Basic checks for required columns
  required_cols <- c("submitter_id", "overall_survival", "deceased")
  if (!all(required_cols %in% colnames(clin))) {
    stop("Clinical file must contain: submitter_id, overall_survival, deceased")
  }
  
  # Clean expression column IDs and match to clinical
  colnames(expr) <- fix_sample_ids(colnames(expr))
  
  # Normalize clinical submitter IDs to 12 chars (safety)
  clin$submitter_id     <- substr(as.character(clin$submitter_id), 1, 12)
  # Ensure numeric survival/time columns
  clin$overall_survival <- as.numeric(clin$overall_survival)
  clin$deceased         <- as.numeric(clin$deceased)
  
  # Match samples between expression and clinical by submitter_id
  shared <- intersect(colnames(expr), clin$submitter_id)
  if (length(shared) < 20) {
    stop(sprintf("Too few matched samples (%d): check sample IDs and inputs.", length(shared)))
  }
  expr <- expr[, shared, drop = FALSE]
  clin <- clin[match(shared, clin$submitter_id), ]
  
  # Return wide patient-level dataframe
  #One row per patient: columns=genes,OS,event
  data.frame(
    t(expr),
    overall_survival = clin$overall_survival,
    deceased         = clin$deceased,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
}

#Apply function
df_LUAD <- prepare_survival_df_precomputed("LUAD_expr_matched.csv", "clinical_LUAD_clean.csv")
df_LUSC <- prepare_survival_df_precomputed("LUSC_expr_matched.csv", "clinical_LUSC_clean.csv")


# Fit a Cox model for one signature (2 or 3 genes) and return summary metrics
evaluate_signature <- function(df, genes) {
  # Ensure genes are in the dataframe
  if (!all(genes %in% colnames(df))) return(NULL)
  subdf <- df[, c(genes, "overall_survival", "deceased"), drop = FALSE]
  subdf <- na.omit(subdf)
  
  # Skip if too few patients, no events, or a gene has zero variance
  if (nrow(subdf) < 20 || length(unique(subdf$deceased)) < 2) return(NULL)
  if (any(sapply(subdf[, genes, drop = FALSE], sd, na.rm = TRUE) == 0)) return(NULL)
  
  # Build Cox formula: Surv ~ geneA + geneB (+ geneC)
  form <- as.formula(
    sprintf("Surv(overall_survival, deceased) ~ %s",
            paste(sprintf("`%s`", genes), collapse = " + "))
  )
  
  fit <- tryCatch(coxph(form, data = subdf), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  #Extract C-index, AIC, logLIK and global wald p-value
  sum_fit <- summary(fit)
  # robust pick of global Wald p-value
  wald_vec <- as.numeric(sum_fit$waldtest)
  wald_p   <- if (!is.null(names(sum_fit$waldtest)) && "pvalue" %in% names(sum_fit$waldtest)) {
    as.numeric(sum_fit$waldtest["pvalue"])
  } else if (length(wald_vec) >= 3) {
    wald_vec[3]
  } else {
    NA_real_
  }
  
  data.frame(
    signature = paste(genes, collapse = " + "),
    C_index   = as.numeric(sum_fit$concordance[1]),
    AIC       = as.numeric(AIC(fit)),
    logLik    = as.numeric(logLik(fit)),
    p_value   = as.numeric(wald_p),
    stringsAsFactors = FALSE
  )
}

# Evaluate all 2- and 3-gene combinations from a candidate list
evaluate_all_combinations <- function(df, top_genes) {
  # Guard against too-short lists
  if (length(top_genes) < 2) stop("Need at least 2 genes to form a signature.")
  pairs   <- combn(top_genes, 2, simplify = FALSE)
  triples <- if (length(top_genes) >= 3) combn(top_genes, 3, simplify = FALSE) else list()
  
  # Fit and collect metrics
  all_sets <- c(pairs, triples)
  res_list <- lapply(all_sets, function(g) evaluate_signature(df, g))
  results  <- do.call(rbind, res_list[!sapply(res_list, is.null)])
  
  if (is.null(results) || nrow(results) == 0) {
    stop("No valid combinations found; adjust filters or sample size.")
  }
  
  #Standardizes numeric types 
  results$C_index <- as.numeric(results$C_index)
  results$AIC     <- as.numeric(results$AIC)
  results$p_value <- as.numeric(results$p_value)
  
  # Rank by discrimination first (C-index desc), then parsimony (AIC asc)
  results[order(-results$C_index, results$AIC), ]
}

# Read Top-15 list and keep only genes present in df
# (strip versions so names line up with df, which uses version-less Ensembl IDs)
top15_LUAD_clean <- gsub("\\..*", "", top15_LUAD)
top15_LUSC_clean <- gsub("\\..*", "", top15_LUSC)
top15_LUAD_in_df <- intersect(top15_LUAD_clean, colnames(df_LUAD))
top15_LUSC_in_df <- intersect(top15_LUSC_clean, colnames(df_LUSC))

# Rank all 2- and 3-gene signatures
ranked_LUAD <- evaluate_all_combinations(df_LUAD, top15_LUAD_in_df)
ranked_LUSC <- evaluate_all_combinations(df_LUSC, top15_LUSC_in_df)

# Compute KM log-rank p-value after turning the signature into a risk score
# Risk score = Cox linear predictor; split by median into High/Low risk groups
get_km_pval <- function(df, signature) {
  genes <- strsplit(signature, " + ", fixed = TRUE)[[1]]
  if (!all(genes %in% colnames(df))) return(1)
  
  subdf <- df[, c(genes, "overall_survival", "deceased"), drop = FALSE]
  subdf <- na.omit(subdf)
  if (nrow(subdf) < 10 || length(unique(subdf$deceased)) < 2) return(1)
  
  form <- as.formula(
    sprintf("Surv(overall_survival, deceased) ~ %s",
            paste(sprintf("`%s`", genes), collapse = " + "))
  )
  fit <- tryCatch(coxph(form, data = subdf), error = function(e) NULL)
  if (is.null(fit)) return(1)
  
  #Splits into High/Low by the median risk score
  risk <- predict(fit, newdata = subdf, type = "lp")
  if (all(is.na(risk)) || sd(risk, na.rm = TRUE) == 0) return(1)
  subdf$risk_group <- factor(risk >= median(risk, na.rm = TRUE), labels = c("Low", "High"))
  if (length(unique(subdf$risk_group)) < 2) return(1)
  
  sdiff <- tryCatch(survdiff(Surv(overall_survival, deceased) ~ risk_group, data = subdf),
                    error = function(e) NULL)
  if (is.null(sdiff)) return(1)
  pval  <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  if (is.na(pval)) pval <- 1
  pval
}

# Compute KM p-values for each signature’s High/Low risk split
ranked_LUAD$km_pval <- sapply(ranked_LUAD$signature, function(sig) get_km_pval(df_LUAD, sig))
ranked_LUSC$km_pval <- sapply(ranked_LUSC$signature, function(sig) get_km_pval(df_LUSC, sig))


# Plot and save KM curves for significant signatures
plot_top_signatures <- function(df, ranked, master_gene_label, prefix, top_n = 5) {
  if (nrow(ranked) == 0) return(invisible(NULL))
  top_n <- min(top_n, nrow(ranked))
  for (i in seq_len(top_n)) {
    genes <- strsplit(ranked$signature[i], " + ", fixed = TRUE)[[1]]
    
    # Map Ensembl IDs to symbols for a nicer title
    gene_syms <- sapply(genes, function(g) {
      sym <- master_gene_label[g]
      if (is.na(sym) || sym == "") g else sym
    })
    
    subdf <- df[, c(genes, "overall_survival", "deceased"), drop = FALSE]
    subdf <- na.omit(subdf)
    if (nrow(subdf) < 10 || length(unique(subdf$deceased)) < 2) next
    
    form <- as.formula(
      sprintf("Surv(overall_survival, deceased) ~ %s",
              paste(sprintf("`%s`", genes), collapse = " + "))
    )
    fit  <- tryCatch(coxph(form, data = subdf), error = function(e) NULL)
    if (is.null(fit)) next
    
    risk <- predict(fit, newdata = subdf, type = "lp")
    if (all(is.na(risk)) || sd(risk, na.rm = TRUE) == 0) next
    subdf$risk_group <- factor(risk >= median(risk, na.rm = TRUE), labels = c("Low", "High"))
    if (length(unique(subdf$risk_group)) < 2) next
    
    km <- survfit(Surv(overall_survival, deceased) ~ risk_group, data = subdf)
    title <- sprintf("%s: %s", prefix, paste(gene_syms, collapse = ", "))
    
    plt <- ggsurvplot(
      km, data = subdf,
      pval       = TRUE,
      risk.table = FALSE,
      palette    = c("#F4A8A8", "#A80000"),
      title      = title,
      ggtheme    = theme_survminer() +
        theme(
          plot.title   = element_text(hjust = 0.5, margin = margin(t = 10, b = 5)),
          plot.margin  = margin(t = 20, r = 20, b = 20, l = 20)
        )
    )
    
    ggsave(
      filename = sprintf("%s_signature_%d_KM.png", tolower(prefix), i),
      plot     = plt$plot,
      width    = 8, height = 5, dpi = 300
    )
  }
}

# Save significant signatures and top 5 by KM p-value
sig_LUAD <- subset(ranked_LUAD, km_pval < 0.05)
sig_LUSC <- subset(ranked_LUSC, km_pval < 0.05)

sig_LUAD <- sig_LUAD[order(sig_LUAD$km_pval), ]
sig_LUSC <- sig_LUSC[order(sig_LUSC$km_pval), ]

#Save results
write.csv(sig_LUAD, "LUAD_top_significant_by_KM_pval.csv", row.names = FALSE)
write.csv(head(sig_LUAD, 5), "LUAD_top5_significant_by_KM_pval.csv", row.names = FALSE)

write.csv(sig_LUSC, "LUSC_top_significant_by_KM_pval.csv", row.names = FALSE)
write.csv(head(sig_LUSC, 5), "LUSC_top5_significant_by_KM_pval.csv", row.names= FALSE)

# Plot significant signatures (KM curves)
if (nrow(sig_LUAD) > 0) {
  plot_top_signatures(df_LUAD, sig_LUAD, master_gene_label, "LUAD")
}
if (nrow(sig_LUSC) > 0) {
  plot_top_signatures(df_LUSC, sig_LUSC, master_gene_label, "LUSC")
}





#=============================mRNA analyses=====================================
###======================== RESEARCH PROJECT: mRNA (Protein-Coding) =============
#---------------- DOWNLOADING DATA FOR LUAD AND LUSC (Tumor) -------------------
# Download tumor data for LUAD and LUSC that will contain the gene expression
# counts per gene (STAR - Counts). We will later subset to protein-coding genes.

# Loading required libraries
library(TCGAbiolinks)      # Querying and downloading TCGA data
library(SummarizedExperiment)

#list and summarize projects
#gdcprojects <- getGDCprojects()
#getProjectSummary("TCGA-LUAD")
#getProjectSummary("TCGA-LUSC")

# Query LUAD Primary Tumor RNA-Seq counts
query_LUAD_mRNA <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# Query LUSC Primary Tumor RNA-Seq counts
query_LUSC_mRNA <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# Download
GDCdownload(query_LUAD_mRNA)
GDCdownload(query_LUSC_mRNA)

# Prepare (SummarizedExperiment: assays, colData, rowData)
LUAD_mRNA_tumor_SE <- GDCprepare(query_LUAD_mRNA)
LUSC_mRNA_tumor_SE <- GDCprepare(query_LUSC_mRNA)


#--------------------------- DOWNLOAD NORMAL SAMPLES ---------------------------
# We also download Normal (Solid Tissue Normal) for both cancer types to allow
# normalization and differential expression between Tumor vs Normal.

query_LUAD_mRNA_ctrl <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

query_LUSC_mRNA_ctrl <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

GDCdownload(query_LUAD_mRNA_ctrl)
GDCdownload(query_LUSC_mRNA_ctrl)

LUAD_mRNA_normal_SE <- GDCprepare(query_LUAD_mRNA_ctrl)
LUSC_mRNA_normal_SE <- GDCprepare(query_LUSC_mRNA_ctrl)


#-------------------- SELECT PROTEIN-CODING GENES (Ensembl) --------------------
# Keep only protein-coding genes (to mirror the feature selection step you did
# for lncRNAs, but now with the protein-coding biotype).

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(dplyr)

# Get all Ensembl genes and subset to protein_coding
all_genes_df <- genes(EnsDb.Hsapiens.v86, return.type = "data.frame")
pc_gene_ids <- all_genes_df %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(gene_id)

# Helper: strip version from Ensembl IDs (e.g., ENSG...1.2 -> ENSG...1)
strip_ver <- function(x) sub("\\..*", "", x)

# Ensure rownames are plain Ensembl IDs (no version)
rownames(LUAD_mRNA_tumor_SE)  <- strip_ver(rownames(LUAD_mRNA_tumor_SE))
rownames(LUSC_mRNA_tumor_SE)  <- strip_ver(rownames(LUSC_mRNA_tumor_SE))
rownames(LUAD_mRNA_normal_SE) <- strip_ver(rownames(LUAD_mRNA_normal_SE))
rownames(LUSC_mRNA_normal_SE) <- strip_ver(rownames(LUSC_mRNA_normal_SE))

# Filter to protein-coding features
LUAD_mRNA_tumor_counts  <- assay(LUAD_mRNA_tumor_SE)[rownames(LUAD_mRNA_tumor_SE) %in% strip_ver(pc_gene_ids), ]
LUSC_mRNA_tumor_counts  <- assay(LUSC_mRNA_tumor_SE)[rownames(LUSC_mRNA_tumor_SE) %in% strip_ver(pc_gene_ids), ]
LUAD_mRNA_normal_counts <- assay(LUAD_mRNA_normal_SE)[rownames(LUAD_mRNA_normal_SE) %in% strip_ver(pc_gene_ids), ]
LUSC_mRNA_normal_counts <- assay(LUSC_mRNA_normal_SE)[rownames(LUSC_mRNA_normal_SE) %in% strip_ver(pc_gene_ids), ]

# Save just the (raw) protein-coding normal counts for downstream matching (optional)
write.csv(LUAD_mRNA_normal_counts, "LUAD_control_mRNA.csv", row.names = TRUE)
write.csv(LUSC_mRNA_normal_counts, "LUSC_control_mRNA.csv", row.names = TRUE)


#-------------------- PREPARE DATASETS FOR NORMALIZATION -----------------------
# Match lncRNA style: fix colnames, add "_tumor"/"_normal", combine.

fix_colnames <- function(df) {
  colnames(df) <- gsub("\\.", "-", colnames(df))
  df
}

make_unique_colnames <- function(df, label) {
  colnames(df) <- paste0(colnames(df), "_", label)
  df
}

LUAD_mRNA_tumor_counts   <- fix_colnames(as.data.frame(LUAD_mRNA_tumor_counts))
LUSC_mRNA_tumor_counts   <- fix_colnames(as.data.frame(LUSC_mRNA_tumor_counts))
LUAD_mRNA_normal_counts  <- fix_colnames(as.data.frame(LUAD_mRNA_normal_counts))
LUSC_mRNA_normal_counts  <- fix_colnames(as.data.frame(LUSC_mRNA_normal_counts))

LUAD_mRNA_tumor_counts   <- make_unique_colnames(LUAD_mRNA_tumor_counts,  "tumor")
LUAD_mRNA_normal_counts  <- make_unique_colnames(LUAD_mRNA_normal_counts, "normal")
LUSC_mRNA_tumor_counts   <- make_unique_colnames(LUSC_mRNA_tumor_counts,  "tumor")
LUSC_mRNA_normal_counts  <- make_unique_colnames(LUSC_mRNA_normal_counts, "normal")

# Counts per group
cat("LUAD Tumor samples:", ncol(LUAD_mRNA_tumor_counts), "\n")
cat("LUAD Normal samples:", ncol(LUAD_mRNA_normal_counts), "\n")
cat("LUSC Tumor samples:", ncol(LUSC_mRNA_tumor_counts), "\n")
cat("LUSC Normal samples:", ncol(LUSC_mRNA_normal_counts), "\n")

# Combine Tumor + Normal (per cancer type)
LUAD_mRNA_combined <- cbind(LUAD_mRNA_tumor_counts, LUAD_mRNA_normal_counts)
LUSC_mRNA_combined <- cbind(LUSC_mRNA_tumor_counts, LUSC_mRNA_normal_counts)

# Save combined raw matrices (used by “Raw vs VST” QC plot)
write.csv(LUAD_mRNA_combined, "LUAD_mRNA_combined.csv")
write.csv(LUSC_mRNA_combined, "LUSC_mRNA_combined.csv")


#============== Check if Normal and Tumor tissues are patient-matched ==========
# Create patient IDs to examine matched pairs (first 3 barcode fields)
get_patient_id <- function(barcode) {
  barcode <- sub("_(tumor|normal)$", "", barcode)
  barcode <- gsub("\\.", "-", barcode)
  parts <- strsplit(barcode, "-")[[1]]
  paste(parts[1:3], collapse = "-")
}

luad_tumor_ids  <- sapply(colnames(LUAD_mRNA_tumor_counts),  get_patient_id)
luad_normal_ids <- sapply(colnames(LUAD_mRNA_normal_counts), get_patient_id)
lusc_tumor_ids  <- sapply(colnames(LUSC_mRNA_tumor_counts),  get_patient_id)
lusc_normal_ids <- sapply(colnames(LUSC_mRNA_normal_counts), get_patient_id)

luad_matches <- data.frame(
  Normal_Sample   = colnames(LUAD_mRNA_normal_counts),
  Patient_ID      = luad_normal_ids,
  Matched_To_Tumor= luad_normal_ids %in% luad_tumor_ids
)
lusc_matches <- data.frame(
  Normal_Sample   = colnames(LUSC_mRNA_normal_counts),
  Patient_ID      = lusc_normal_ids,
  Matched_To_Tumor= lusc_normal_ids %in% lusc_tumor_ids
)

cat("LUAD matched normals to tumors:\n")
print(table(luad_matches$Matched_To_Tumor))
cat("LUSC matched normals to tumors:\n")
print(table(lusc_matches$Matched_To_Tumor))


#------------------------------- Normalize counts ------------------------------
# Run DESeq2 with design ~ condition (Tumor vs Normal), keep genes with at least
# 10 counts in 5 or more samples (to mirror lncRNA filtering).

library(DESeq2)

LUAD_condition <- factor(c(rep("Tumor", ncol(LUAD_mRNA_tumor_counts)),
                           rep("Normal", ncol(LUAD_mRNA_normal_counts))))
LUSC_condition <- factor(c(rep("Tumor", ncol(LUSC_mRNA_tumor_counts)),
                           rep("Normal", ncol(LUSC_mRNA_normal_counts))))

LUAD_colData <- data.frame(condition = LUAD_condition, row.names = colnames(LUAD_mRNA_combined))
LUSC_colData <- data.frame(condition = LUSC_condition, row.names = colnames(LUSC_mRNA_combined))

dds_LUAD_mRNA <- DESeqDataSetFromMatrix(countData = LUAD_mRNA_combined,
                                        colData = LUAD_colData, design = ~ condition)
dds_LUSC_mRNA <- DESeqDataSetFromMatrix(countData = LUSC_mRNA_combined,
                                        colData = LUSC_colData, design = ~ condition)

keep_LUAD <- rowSums(counts(dds_LUAD_mRNA) >= 10) >= 5
keep_LUSC <- rowSums(counts(dds_LUSC_mRNA) >= 10) >= 5
dds_LUAD_mRNA <- dds_LUAD_mRNA[keep_LUAD, ]
dds_LUSC_mRNA <- dds_LUSC_mRNA[keep_LUSC, ]

dds_LUAD_mRNA <- DESeq(dds_LUAD_mRNA)
dds_LUSC_mRNA <- DESeq(dds_LUSC_mRNA)

# Save DESeq2 objects if needed
saveRDS(dds_LUAD_mRNA, "dds_LUAD_mRNA.rds")
saveRDS(dds_LUSC_mRNA, "dds_LUSC_mRNA.rds")

# Get Normalized counts
LUAD_mRNA_norm_counts <- counts(dds_LUAD_mRNA, normalized = TRUE)
LUSC_mRNA_norm_counts <- counts(dds_LUSC_mRNA, normalized = TRUE)
write.csv(LUAD_mRNA_norm_counts, "LUAD_mRNA_normalized_counts_DESeq2.csv")
write.csv(LUSC_mRNA_norm_counts, "LUSC_mRNA_normalized_counts_DESeq2.csv")

# Log-transform (VST) for visualization
LUAD_mRNA_vst <- vst(dds_LUAD_mRNA, blind = TRUE)
LUSC_mRNA_vst <- vst(dds_LUSC_mRNA, blind = TRUE)

# Preserve the tumor/normal suffix labels
vst_LUAD <- assay(LUAD_mRNA_vst); colnames(vst_LUAD) <- colnames(LUAD_mRNA_combined)
vst_LUSC <- assay(LUSC_mRNA_vst); colnames(vst_LUSC) <- colnames(LUSC_mRNA_combined)

write.csv(vst_LUAD, "LUAD_mrna_log_trans_after_norm.csv")
write.csv(vst_LUSC, "LUSC_mrna_log_trans_after_norm.csv")

# Also save VST of NORMAL-ONLY (so later clinical matching for normals can load it)
LUAD_vst_normal_only <- vst_LUAD[, grepl("_normal$", colnames(vst_LUAD), ignore.case = TRUE), drop = FALSE]
LUSC_vst_normal_only <- vst_LUSC[, grepl("_normal$", colnames(vst_LUSC), ignore.case = TRUE), drop = FALSE]
write.csv(LUAD_vst_normal_only, "LUAD_control_mRNA_log_trans_after_norm.csv")
write.csv(LUSC_vst_normal_only, "LUSC_control_mRNA_log_trans_after_norm.csv")

#----------------------------- Differential Expression -------------------------
# Perform DEA Tumor vs Normal. Save both unshrunken and shrunken results,
# and also a filtered shrunken table (padj < 0.01 & |log2FC| > 1).

res_LUAD_mRNA <- results(dds_LUAD_mRNA, contrast = c("condition", "Tumor", "Normal"))
res_LUSC_mRNA <- results(dds_LUSC_mRNA, contrast = c("condition", "Tumor", "Normal"))

library(apeglm)
res_LUAD_mRNA_shrink <- lfcShrink(dds_LUAD_mRNA, coef = "condition_Tumor_vs_Normal", type = "apeglm")
res_LUSC_mRNA_shrink <- lfcShrink(dds_LUSC_mRNA, coef = "condition_Tumor_vs_Normal", type = "apeglm")

# Filtered shrunken (strict cutoff as used in plots)
res_LUAD_mRNA_filtered <- subset(as.data.frame(res_LUAD_mRNA_shrink), padj < 0.01 & abs(log2FoldChange) > 1)
res_LUSC_mRNA_filtered <- subset(as.data.frame(res_LUSC_mRNA_shrink), padj < 0.01 & abs(log2FoldChange) > 1)

# Save full + shrunken + filtered
write.csv(as.data.frame(res_LUAD_mRNA),        "LUAD_mrna_DESeq2_all_results.csv")
write.csv(as.data.frame(res_LUAD_mRNA_shrink), "LUAD_mrna_DESeq2_all_results_shrunk.csv")
write.csv(res_LUAD_mRNA_filtered,              "LUAD_mrna_DESeq2_filtered_shrunk.csv")

write.csv(as.data.frame(res_LUSC_mRNA),        "LUSC_mrna_DESeq2_all_results.csv")
write.csv(as.data.frame(res_LUSC_mRNA_shrink), "LUSC_mrna_DESeq2_all_results_shrunk.csv")
write.csv(res_LUSC_mRNA_filtered,              "LUSC_mrna_DESeq2_filtered_shrunk.csv")


#---------------------------- MA plots to visualize DEA ------------------------
library(ggplot2)

save_ma_plot <- function(res_df, cancer_type, out_png) {
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df$significant <- res_df$padj < 0.05
  
  p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "#F4A8A8", "TRUE" = "#A80000")) +
    scale_x_log10() +
    labs(
      title = paste("MA Plot -", cancer_type, "(mRNA)"),
      x = "Mean Expression (log10)",
      y = "log2 Fold Change",
      color = "padj < 0.05"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 28),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text    = element_text(size = 24),
      legend.title = element_text(size = 24),
      legend.text  = element_text(size = 24),
      legend.position = "bottom"
    )
  ggsave(out_png, plot = p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

save_ma_plot(as.data.frame(res_LUAD_mRNA_shrink), "LUAD", "MA_plot_LUAD_mRNA.png")
save_ma_plot(as.data.frame(res_LUSC_mRNA_shrink), "LUSC", "MA_plot_LUSC_mRNA.png")


#---------------- Map Ensembl IDs to HGNC Symbols using biomaRt ----------------
library(biomaRt)

deg_LUAD_mRNA_allshrink <- read.csv("LUAD_mrna_DESeq2_all_results_shrunk.csv", row.names = 1)
deg_LUSC_mRNA_allshrink <- read.csv("LUSC_mrna_DESeq2_all_results_shrunk.csv", row.names = 1)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_ensembl_ids <- unique(c(rownames(deg_LUAD_mRNA_allshrink), rownames(deg_LUSC_mRNA_allshrink)))
all_ensembl_ids <- strip_ver(all_ensembl_ids)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_ensembl_ids,
  mart = mart
)

write.csv(gene_map, "gene_map_mRNA_LUAD_LUSC.csv", row.names = FALSE)
master_gene_label <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)
save(master_gene_label, file = "master_gene_label_mRNA.RData")


#---------------------------- Select Top 15 DE mRNAs ---------------------------
library(dplyr)

sig_deg_LUAD_mRNA <- deg_LUAD_mRNA_allshrink %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange)))

sig_deg_LUSC_mRNA <- deg_LUSC_mRNA_allshrink %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange)))

top15_LUAD_mRNA <- rownames(sig_deg_LUAD_mRNA)[1:min(15, nrow(sig_deg_LUAD_mRNA))]
top15_LUSC_mRNA <- rownames(sig_deg_LUSC_mRNA)[1:min(15, nrow(sig_deg_LUSC_mRNA))]

write.csv(data.frame(gene_id = top15_LUAD_mRNA), "LUAD_top15_mRNA_DEGs.csv", row.names = FALSE)
write.csv(data.frame(gene_id = top15_LUSC_mRNA), "LUSC_top15_mRNA_DEGs.csv", row.names = FALSE)

cat("Top 15 LUAD mRNA genes:\n"); print(top15_LUAD_mRNA)
cat("\nTop 15 LUSC mRNA genes:\n"); print(top15_LUSC_mRNA)

# Also save “with gene names” tables (for your correlation section)
LUAD_mRNA_label_map <- subset(gene_map, ensembl_gene_id %in% strip_ver(top15_LUAD_mRNA))
LUSC_mRNA_label_map <- subset(gene_map, ensembl_gene_id %in% strip_ver(top15_LUSC_mRNA))
LUAD_mRNA_label_map$hgnc_symbol[is.na(LUAD_mRNA_label_map$hgnc_symbol) | LUAD_mRNA_label_map$hgnc_symbol == ""] <- LUAD_mRNA_label_map$ensembl_gene_id[is.na(LUAD_mRNA_label_map$hgnc_symbol) | LUAD_mRNA_label_map$hgnc_symbol == ""]
LUSC_mRNA_label_map$hgnc_symbol[is.na(LUSC_mRNA_label_map$hgnc_symbol) | LUSC_mRNA_label_map$hgnc_symbol == ""] <- LUSC_mRNA_label_map$ensembl_gene_id[is.na(LUSC_mRNA_label_map$hgnc_symbol) | LUSC_mRNA_label_map$hgnc_symbol == ""]
write.csv(LUAD_mRNA_label_map, "LUAD_top15_mRNA_with_gene_names.csv", row.names = FALSE)
write.csv(LUSC_mRNA_label_map, "LUSC_top15_mRNA_with_gene_names.csv", row.names = FALSE)


#------------------------- Boxplots of Top 15 DE mRNAs -------------------------
library(ggpubr)
library(reshape2)

# Load VST (mRNA) and replace rownames column
LUAD_vst_mRNA <- read.csv("LUAD_mrna_log_trans_after_norm.csv", row.names = 1)
LUSC_vst_mRNA <- read.csv("LUSC_mrna_log_trans_after_norm.csv", row.names = 1)

LUAD_vst_mRNA$Gene <- strip_ver(rownames(LUAD_vst_mRNA))
LUSC_vst_mRNA$Gene <- strip_ver(rownames(LUSC_vst_mRNA))

# Subset top 15
LUAD_top_mRNA <- LUAD_vst_mRNA[LUAD_vst_mRNA$Gene %in% strip_ver(top15_LUAD_mRNA), ]
LUSC_top_mRNA <- LUSC_vst_mRNA[LUSC_vst_mRNA$Gene %in% strip_ver(top15_LUSC_mRNA), ]

# Melt
LUAD_melt_mRNA <- melt(LUAD_top_mRNA, id.vars = "Gene")
LUSC_melt_mRNA <- melt(LUSC_top_mRNA, id.vars = "Gene")

# SampleType
LUAD_melt_mRNA$SampleType <- ifelse(grepl("_tumor$", LUAD_melt_mRNA$variable, ignore.case = TRUE), "Tumor", "Normal")
LUSC_melt_mRNA$SampleType <- ifelse(grepl("_tumor$", LUSC_melt_mRNA$variable, ignore.case = TRUE), "Tumor", "Normal")

# Replace Ensembl IDs with symbols when available
LUAD_melt_mRNA$Gene <- ifelse(
  is.na(master_gene_label[LUAD_melt_mRNA$Gene]) | master_gene_label[LUAD_melt_mRNA$Gene] == "",
  LUAD_melt_mRNA$Gene,
  master_gene_label[LUAD_melt_mRNA$Gene]
)
LUSC_melt_mRNA$Gene <- ifelse(
  is.na(master_gene_label[LUSC_melt_mRNA$Gene]) | master_gene_label[LUSC_melt_mRNA$Gene] == "",
  LUSC_melt_mRNA$Gene,
  master_gene_label[LUSC_melt_mRNA$Gene]
)

# Plot helper
plot_top15_boxplot <- function(data, title, filename) {
  p <- ggplot(data, aes(x = Gene, y = value, fill = SampleType)) +
    geom_boxplot(width = 0.85, outlier.shape = NA, alpha = 0.7) +
    stat_compare_means(aes(group = SampleType),
                       method = "wilcox.test",
                       label = "p.signif",
                       size = 12,
                       tip.length = 0.01,
                       label.y = max(data$value) * 1.2) +
    scale_fill_manual(values = c("Normal" = "#F4A8A8", "Tumor" = "#A80000")) +
    coord_cartesian(ylim = c(min(data$value), max(data$value) * 1.3)) +
    theme_minimal(base_size = 18) +
    ggtitle(title) +
    ylab("VST") + xlab("Gene") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y  = element_text(size = 20),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text  = element_text(size = 20),
      legend.position = "bottom"
    )
  ggsave(filename, plot = p, width = 16, height = 10, dpi = 300)
  p
}

p_luad_mrna <- plot_top15_boxplot(LUAD_melt_mRNA, "Top 15 Differentially Expressed mRNAs in LUAD", "LUAD_Top15_mRNA_Boxplot.png")
p_lusc_mrna <- plot_top15_boxplot(LUSC_melt_mRNA, "Top 15 Differentially Expressed mRNAs in LUSC", "LUSC_Top15_mRNA_Boxplot.png")

print(p_luad_mrna)
print(p_lusc_mrna)


#----------------------- Download clinical data (same style) -------------------
library(readr)

clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD")
clinical_LUSC <- GDCquery_clinic(project = "TCGA-LUSC")

# Binary vital + overall survival (days)
clinical_LUAD$deceased <- ifelse(clinical_LUAD$vital_status == "Alive", FALSE, TRUE)
clinical_LUSC$deceased <- ifelse(clinical_LUSC$vital_status == "Alive", FALSE, TRUE)

clinical_LUAD$overall_survival <- ifelse(
  clinical_LUAD$deceased, clinical_LUAD$days_to_death, clinical_LUAD$days_to_last_follow_up
)
clinical_LUSC$overall_survival <- ifelse(
  clinical_LUSC$deceased, clinical_LUSC$days_to_death, clinical_LUSC$days_to_last_follow_up
)

write_csv(clinical_LUAD, "clinical_LUAD_clean.csv")
write_csv(clinical_LUSC, "clinical_LUSC_clean.csv")



#----------------------------- Correlation of mRNAs and lncRNAs-----------------
library(corrplot)

#safe TCGA ID cleaner (keeps TCGA-XX-XXXX if present; never returns NA)
clean_tcga_3_safe <- function(x) {
  x0 <- as.character(x)
  x  <- gsub("\\s+", "", x0)
  x  <- gsub("\\.", "-", x)
  x  <- gsub("(_tumor|_normal)$", "", x, ignore.case = TRUE)
  x  <- toupper(x)
  m    <- regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x)
  keep <- ifelse(m > 0, regmatches(x, m), NA_character_)
  fb <- sapply(strsplit(x, "-"), function(p)
    if (length(p) >= 3) paste(p[1:3], collapse = "-") else NA_character_)
  out <- ifelse(!is.na(keep), keep, fb)
  out[is.na(out)] <- x0[is.na(out)]
  out
}

# collapse duplicate columns created by cleaning (mean across dupes)
collapse_cols_by_id <- function(mat) {
  ids <- colnames(mat)
  idx <- split(seq_along(ids), ids)
  if (any(lengths(idx) > 1)) {
    mat <- sapply(idx, function(j) rowMeans(as.matrix(mat[, j, drop = FALSE])))
    mat <- as.matrix(mat)
    colnames(mat) <- names(idx)
  }
  mat
}

#compute correlation and p-value matrices (lncRNA rows x mRNA cols)
cor_and_p <- function(A, B, method = "pearson") {
  stopifnot(ncol(A) == ncol(B))
  # r
  r <- cor(t(A), t(B), method = method, use = "pairwise.complete.obs")
  # p
  p <- matrix(NA_real_, nrow = nrow(A), ncol = nrow(B))
  rownames(p) <- rownames(A); colnames(p) <- rownames(B)
  for (i in seq_len(nrow(A))) {
    ai <- as.numeric(A[i, ])
    for (j in seq_len(nrow(B))) {
      bj <- as.numeric(B[j, ])
      tt <- suppressWarnings(cor.test(ai, bj, method = method))
      p[i, j] <- tt$p.value
    }
  }
  list(r = r, p = p)
}

#run once per cancer
run_one_cancer <- function(
    cancer, 
    lnc_file, mrna_file, 
    top15_lnc_csv, top15_mrna_csv
) {
  message("=== ", cancer, " ===")
  # read expression (keep column names as-is)
  lnc <- read.csv(lnc_file,  row.names = 1, check.names = FALSE)
  mrn <- read.csv(mrna_file, row.names = 1, check.names = FALSE)
  
  #keep tumor-only columns by suffix
  lnc <- lnc[, grepl("_tumor$", colnames(lnc)), drop = FALSE]
  mrn <- mrn[, grepl("_tumor$", colnames(mrn)), drop = FALSE]
  
  # clean sample IDs + collapse duplicates
  colnames(lnc) <- clean_tcga_3_safe(colnames(lnc))
  colnames(mrn) <- clean_tcga_3_safe(colnames(mrn))
  lnc <- collapse_cols_by_id(lnc)
  mrn <- collapse_cols_by_id(mrn)
  
  # align samples
  common <- intersect(colnames(lnc), colnames(mrn))
  message("Matched tumor samples: ", length(common))
  lnc <- lnc[, common, drop = FALSE]
  mrn <- mrn[, common, drop = FALSE]
  
  # read Top15 lists (Ensembl IDs + symbols)
  tlnc <- read.csv(top15_lnc_csv)
  tmrn <- read.csv(top15_mrna_csv)
  lnc_ids <- gsub("\\..*", "", tlnc$ensembl_gene_id)
  mrn_ids <- gsub("\\..*", "", tmrn$ensembl_gene_id)
  lnc_syms <- ifelse(is.na(tlnc$hgnc_symbol) | tlnc$hgnc_symbol == "", lnc_ids, tlnc$hgnc_symbol)
  mrn_syms <- ifelse(is.na(tmrn$hgnc_symbol) | tmrn$hgnc_symbol == "", mrn_ids, tmrn$hgnc_symbol)
  
  # clean rownames in expression, subset to Top15
  rownames(lnc) <- gsub("\\..*", "", rownames(lnc))
  rownames(mrn) <- gsub("\\..*", "", rownames(mrn))
  lnc <- lnc[intersect(lnc_ids, rownames(lnc)), , drop = FALSE]
  mrn <- mrn[intersect(mrn_ids, rownames(mrn)), , drop = FALSE]
  # keep Top15 order
  lnc <- lnc[lnc_ids[lnc_ids %in% rownames(lnc)], , drop = FALSE]
  mrn <- mrn[mrn_ids[mrn_ids %in% rownames(mrn)], , drop = FALSE]
  # set pretty names
  rownames(lnc) <- lnc_syms[match(rownames(lnc), lnc_ids)]
  rownames(mrn) <- mrn_syms[match(rownames(mrn), mrn_ids)]
  
  # correlation + p-values (+ FDR)
  cp <- cor_and_p(lnc, mrn, method = "pearson")
  rmat <- cp$r
  pmat <- cp$p
  padj <- matrix(p.adjust(as.vector(pmat), method = "BH"),
                 nrow = nrow(pmat), ncol = ncol(pmat),
                 dimnames = dimnames(pmat))
  
  #  save CSVs
  write.csv(rmat,  paste0(cancer, "_lncRNA_mRNA_cor.csv"))
  write.csv(pmat,  paste0(cancer, "_lncRNA_mRNA_pval.csv"))
  write.csv(padj,  paste0(cancer, "_lncRNA_mRNA_padj.csv"))
  
  # corrplot (mask by FDR<0.05; non-sig blank)
  png(paste0(cancer, "_lncRNA_mRNA_corrplot.png"), width = 1400, height = 1050, res = 160)
  corrplot(
    rmat,
    method    = "circle",
    type      = "full",
    order     = "hclust",
    tl.col    = "black",
    tl.srt    = 45,
    p.mat     = padj,         # FDR matrix
    sig.level = c(0.001, 0.01, 0.05),
    insig     = "label_sig",  # draw stars for significance
    pch.cex   = 1.2,          # star size
    pch.col   = "black"       # star color
  )
  
  
  dev.off()
  
  message("Saved files for ", cancer)
}

#Run LUAD
run_one_cancer(
  cancer        = "LUAD",
  lnc_file      = "LUAD_lncRNA_log_trans_after_norm.csv",
  mrna_file     = "LUAD_mrna_log_trans_after_norm.csv",
  top15_lnc_csv = "LUAD_top15_lncRNA_with_gene_names.csv",
  top15_mrna_csv= "LUAD_top15_mRNA_with_gene_names.csv"
)

#Run LUSC
run_one_cancer(
  cancer        = "LUSC",
  lnc_file      = "LUSC_lncRNA_log_trans_after_norm.csv",
  mrna_file     = "LUSC_mrna_log_trans_after_norm.csv",
  top15_lnc_csv = "LUSC_top15_lncRNA_with_gene_names.csv",
  top15_mrna_csv= "LUSC_top15_mRNA_with_gene_names.csv"
)




