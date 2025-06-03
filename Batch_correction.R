BiocManager::install("sva")
install.packages("devtools")
devtools::install_github("zhangyuqing/sva-devel")
library(sva)
library(edgeR)
library(ggplot2)
library(ggfortify) 
library(gridExtra)

###Data###
tcga_counts <- read.csv("Common-count-TCGA.csv", row.names = 1)
external_counts <- read.csv("Common-count-SRA.csv", row.names = 1)
merged_counts <- cbind(tcga_counts, external_counts)
merged_counts <- cbind(tcga_counts, external_counts)
head(merged_counts, 5)
dge <- DGEList(counts = merged_counts)
dge <- calcNormFactors(dge, method = "TMM")
expression_data <- log2(cpm(dge) + 1)
keep_genes <- apply(expression_data, 1, var) > 0
expression_data <- expression_data[keep_genes, ]

######### Group info ##########
group_info <- read.csv("Group.csv")
group_info <- group_info[match(colnames(expression_data), group_info$Sample), ]

# Create model matrix for biological variable of interest (Tissue)
mod <- model.matrix(~Tissue, data = group_info)

# Batch correction with ComBat
# Using both batch and preserving biological difference (Tissue)
batch_corrected_data <- ComBat(dat = as.matrix(expression_data),
                               batch = group_info$Batch,
                               mod = mod,
                               par.prior = TRUE,
                               prior.plots = FALSE)
# Convert back to data frame
batch_corrected_data <- as.data.frame(batch_corrected_data)
# Save the data
write.csv(batch_corrected_data, "batch_corrected_data_1.csv")

################ Plot #####################
library(ggplot2)
library(gridExtra)

plot_pca <- function(data, title, group_info) {
  pca_res <- prcomp(t(data))
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Batch <- factor(group_info$Batch)
  pca_df$Tissue <- factor(group_info$Tissue)
  
  # Create a combined Batch-Tissue variable for styling
  pca_df$Batch_Tissue <- interaction(pca_df$Batch, pca_df$Tissue)
  
  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(shape = Batch, color = Batch_Tissue), size = 2) +
    ggtitle(title) +
    theme_minimal() +
    # Shapes: Different shape for each batch
    scale_shape_manual(values = c(16, 17), name = "Batch") + # 16=circle, 17=triangle
    # Colors: Different colors for tumor/normal within each batch
    scale_color_manual(
      values = c(
        "TCGA-STAD.Tumor" = "cadetblue1",  # Red for Batch1-Tumor
        "TCGA-STAD.Normal" = "green3",  # Green for Batch1-Normal
        "PRJNA435914.Tumor" = "brown2",   # Blue for Batch2-Tumor
        "PRJNA435914.Normal" = "yellow3"   # Purple for Batch2-Normal
      ),
      name = "Tissue",
      labels = c("TCGA-STAD Normal", "PRJNA435914 Normal", "PRJNA435914 Tumor", "TCGA-STAD Tumor"),
      guide = guide_legend(override.aes = list(shape = 16))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    labs(
      x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(color = "black")),
      color = guide_legend(order = 2)
    )
}

# Generate plots
p_before <- plot_pca(expression_data, "PCA Before Batch Correction", group_info)
p_after <- plot_pca(batch_corrected_data, "PCA After Batch Correction", group_info)

# Display plots
grid.arrange(p_before, p_after, ncol = 2)


############################ Two color ###############################
plot_pca <- function(data, title, group_info) {
  pca_res <- prcomp(t(data))
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Batch <- factor(group_info$Batch)
  pca_df$Tissue <- factor(group_info$Tissue)
  
  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(shape = Batch, color = Tissue), size = 2) +  # Color by Tissue only
    ggtitle(title) +
    theme_minimal() +
    # Shapes: Different shape for each batch (circle for TCGA, triangle for PRJNA)
    scale_shape_manual(values = c(16, 1), name = "Batch") +
    # Colors: Only two colors for Tumor/Normal (ignore batch in coloring)
    scale_color_manual(
      values = c(
        "Tumor" = "red1",      # Same red for all tumor samples
        "Normal" = "blue4"    # Same green for all normal samples
      ),
      name = "Tissue"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    labs(
      x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(color = "black")),
      color = guide_legend(order = 2)
    )  # This was the missing closing parenthesis
}

# Generate plots
p_before <- plot_pca(expression_data, "PCA Before Batch Correction", group_info)
p_after <- plot_pca(batch_corrected_data, "PCA After Batch Correction", group_info)

# Display plots
grid.arrange(p_before, p_after, ncol = 2)
