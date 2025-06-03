# Install Bioconductor package manager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install DESeq2 if not installed:
BiocManager::install("DESeq2")

# Load required libraries:
library(DESeq2)

#load the count matrix:
counts <- read.csv("C:/Users/ACM/OneDrive/Desktop/Project/SRA/lncRNA_counts.csv", row.names = 1)  # Set first column as row names (gene IDs)
View(counts)
#load metadata:
metadata <- read.csv("C:/Users/ACM/OneDrive/Desktop/Project/SRA/meta_sra.csv", row.names = 1)  # Set first column as row names (sample IDs)
View(metadata)

#DESeq dataset creation:
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~ Sample)

#Remove Low-Expressed Genes (Filtering)
dds <- dds[rowSums(counts(dds)) > 10, ]

#Normalize the Data:
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

#Perform Differential Expression Analysis:
dds <- DESeq(dds)  # Run the DESeq2 pipeline
res <- results(dds)  # Get the results

#convert to CSV:
write.csv(res, "C:/Users/ACM/OneDrive/Desktop/Project/SRA/DEG Analysis_SRA.csv")



