# MSc Thesis Project: Prognostic lncRNA Biomarker Discovery in Gastric Adenocarcinoma

This repository contains the workflow and results of my MSc research project, which focused on identifying **prognostic long non-coding RNA (lncRNA) biomarkers** in **gastric adenocarcinoma** using RNA-Seq datasets.

## 🎯 Objective
To identify lncRNAs significantly associated with patient survival outcomes using survival analysis on transcriptomic data from gastric cancer patients.

## 📂 Datasets Used
- **TCGA-STAD** (The Cancer Genome Atlas)
- **PRJNA435914** (European Nucleotide Archive)
## 📥 Data Access
Due to large file sizes, raw RNA-Seq data is **not included** in this repository.

To download the data:
- **TCGA-STAD**: Use [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) in R
- **PRJNA435914**: Download from [ENA](https://www.ebi.ac.uk/ena/browser/text-search?query=PRJNA435914) or using `fastq-dump` via SRA Toolkit

## 🔬 Methods
- RNA-Seq data preprocessing (QC, alignment using HISAT2, transcript assembly with StringTie)
- Differential expression analysis using **DESeq2**
- Feature selection using **Elastic Net** and **Boruta**
- **Survival analysis** using Kaplan-Meier estimates and Cox Proportional Hazards models

## 🧠 Tools & Technologies
- **Languages**: Python, R  
- **Bioinformatics Tools**: HISAT2, SAMtools, StringTie  
- **Statistical Packages**: DESeq2, Survival, glmnet, Boruta

## 📈 Key Findings
- Identified **145 common differentially expressed lncRNAs**
- Reduced to **37 high-confidence lncRNAs** via Elastic Net and Boruta
- **3 lncRNAs (LINC01614, PGM5-AS1, FRMD6-AS2)** showed significant association with overall survival  
- These lncRNAs are promising candidates for **prognostic biomarkers** in gastric cancer


## 👩‍💻 Author
**Fahiz Mohammed PP**  
MSc in Bioinformatics  
Pondicherry University

