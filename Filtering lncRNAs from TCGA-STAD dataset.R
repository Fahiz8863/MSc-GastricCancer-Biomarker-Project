# ---------------------------
# Step 1: Load lncRNA Annotation File (Ensembl Gene IDs only)
# ---------------------------

lnc_annot <- read.csv("lncRNA_annotation.csv", stringsAsFactors = FALSE)

# Remove version numbers (e.g., ENSG00000123456.1 → ENSG00000123456)
lnc_annot$gene_id <- sub("\\..*", "", lnc_annot$gene_id)

# Vector of unique Ensembl IDs for lncRNAs
lnc_ids <- unique(lnc_annot$gene_id)

# ---------------------------
# Step 2: Prepare TCGA Expression Matrix from SummarizedExperiment
# ---------------------------

library(SummarizedExperiment)

# Remove version numbers from TCGA rownames (if present)
rownames(expdat) <- sub("\\..*", "", rownames(expdat))

# Extract the expression matrix
expr_matrix <- assay(expdat)

# ---------------------------
# Step 3: Filter for lncRNAs
# ---------------------------

# Keep only rows whose Ensembl IDs match lncRNA IDs
lnc_expr <- expr_matrix[rownames(expr_matrix) %in% lnc_ids, ]

# ---------------------------
# Step 4: Save Filtered Matrix
# ---------------------------

write.csv(lnc_expr, "TCGA_STAD_lncRNA_expression.csv")
