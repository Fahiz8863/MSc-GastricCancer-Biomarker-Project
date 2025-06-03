# ---------------------------
# Script: annotate_lncRNA_counts.R
# Purpose: Merge lncRNA count matrix with gene names using Ensembl IDs
# Author: [Your Name]
# ---------------------------

# Load required library
library(dplyr)

# ---------------------------
# Step 1: Load Data
# ---------------------------

# Read lncRNA count matrix (with Ensembl IDs)
count_data <- read.csv("C:/Users/HP/Desktop/Fahiz/Project/TCGA/data/lnc Result/LNC_UpandDown.csv", 
                       header = TRUE, stringsAsFactors = FALSE)

# Remove version numbers from Ensembl IDs
count_data$gene_id <- gsub("\\..*", "", count_data$gene_id)

# Read gene annotation file (with Ensembl IDs and gene names)
annotation_data <- read.csv("C:/Users/HP/Desktop/Fahiz/Project/TCGA/data/lnc Result/1_Gene_name.csv", 
                            header = TRUE, stringsAsFactors = FALSE)

# Remove version numbers
annotation_data$gene_id <- gsub("\\..*", "", annotation_data$gene_id)

# ---------------------------
# Step 2: Merge by Ensembl Gene ID
# ---------------------------

merged_data <- merge(count_data, annotation_data, by = "gene_id")

# Reorder columns: gene_id, gene_name, followed by expression data
merged_data <- merged_data[, c("gene_id", "gene_name", colnames(count_data)[-1])]

# Remove gene_id column if not needed
merged_data <- merged_data[, -which(colnames(merged_data) == "gene_id")]

# ---------------------------
# Step 3: Output
# ---------------------------

# Preview the result
head(merged_data)

# Save to CSV
write.csv(merged_data, 
          "C:/Users/HP/Desktop/Fahiz/Project/TCGA/data/lnc Result/DEG_analysis_TCGA.csv", 
          row.names = FALSE)
