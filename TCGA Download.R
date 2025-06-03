``` R

### TCGA Data Downloads ##

## Installation of BioCmanager and TCGAbiolinks ##
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

## Updated Version ##
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

##if not use devtools or remotes##

### Importing Libraries ### 
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(tidyr)

## Project ##
## TCGA-STAD ##

### Query for mRNA ##

## The following options are used to search mRNA results using TCGAbiolinks:
## Harmonized database: data aligned against hg38 Reference Genome
## data.category: "Transcriptome Profiling"
## data.type: "Gene Expression Quantification"
## workflow.type = "STAR - Counts"

query.exp.hg38 <- GDCquery(
  project = "TCGA-STAD", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
)
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(
  query = query.exp.hg38,
  save = TRUE, 
  save.filename = "exp.rda"
)

# Query clinical data
query.clin <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Clinical",
  file.type = "xml"
)

# Download clinical data
GDCdownload(query.clin)

# Prepare the clinical data into a readable format
clinical.data <- GDCprepare_clinic(query.clin, clinical.info = "patient")
