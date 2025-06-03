

# Load necessary libraries
library(ggplot2)

# Read your DEA results file
# Assuming your file has columns: gene_id, log2FoldChange, padj
data <- read.csv("DEG Analysis_SRA.csv", header = TRUE)

# Add a new column for significance status
data$significance <- "Not Significant"
data$significance[data$log2FoldChange > 1 & data$padj < 0.05] <- "Upregulated"
data$significance[data$log2FoldChange < -1 & data$padj < 0.05] <- "Downregulated"

# Plotting the Volcano plot
ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "grey")) +
  labs(title = "Volcano Plot of Differentially Expressed lncRNAs (SRA Data)",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
