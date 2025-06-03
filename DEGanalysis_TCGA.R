library("DESeq2")
library(dplyr)
library('biomaRt')
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(forcats)
library("ggrepel")

## Uploading Data ##

data <- read.csv("C:/Users/ACM/OneDrive/Desktop/Project/TCGA/data/lnc_F.csv", header = T)
meta <- read.csv("C:/Users/ACM/OneDrive/Desktop/Project/TCGA/data/col_meta_F.csv", row.names = 1)

#convert the meta to those with required column:
meta <- meta[, c("barcode", "tissue_type")]

#to convert '-' to '.' :
meta$barcode <- gsub("-", ".", meta$barcode)

#change the raw name of the data:
rownames(data) <- data$Column1
data$Column1<- NULL

#change the raw name of the meta:
rownames(meta) <- meta$barcode
meta$barcode <- NULL 

#checking column names consistency:
all(colnames(data) %in% rownames(meta))

#ensuring the exact order match
all(colnames(data) == rownames(meta))

#loading DESeq2 library:
library("DESeq2")

#creating DESeq2 data set:
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta,
                              design = ~ tissue_type)

#view DESeq dataset:
dds

#filtering low count genes:
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#setting factor levels:
dds$tissue_type <- factor(dds$tissue_type, levels = c("Tumor","Normal"))

#changing the reference level:
dds$tissue_type <- relevel(dds$tissue_type, ref = "Normal")

#dropping unused factor levels:
dds$tissue_type <- droplevels(dds$tissue_type)

#running the DESeq2 pipeline:
dds <- DESeq(dds)

#Extracting differential expression results:
res <- results(dds)
res
summary(res)

write.csv(subset(res), "res.csv")

