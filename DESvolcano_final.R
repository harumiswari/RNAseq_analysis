library(DESeq2)
library(readxl)
library(ggplot2)

# Close any open graphical devices
while(dev.cur() > 1) dev.off()

# Read raw counts data
# raw_counts <- read_excel("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/Raw_Read_Mouse.xlsx", sheet = 1)
raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/rna_readcount_mouse_invivo.csv", header = TRUE, row.names = 1)
raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/Raw_Read_Mouse.csv", header = TRUE, row.names = 1)

# Check for and handle NAs
if (any(is.na(raw_counts))) {
  # If there are NAs, remove rows with NAs
  raw_counts <- raw_counts[complete.cases(raw_counts), ]
}

# Filter out genes with low read counts
filtered_counts <- raw_counts[rowSums(raw_counts) > 10, ]

# Create colData with treatment information
colData <- data.frame(
  sample = colnames(filtered_counts),
  treatment = factor(c(rep("3wk_Ckm", 4), rep("3wk_SCR", 4), rep("8wk_Ckm", 4), rep("8wk_SCR", 3), rep("NT", 2))),
  replicate = factor(c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3","4", "1", "2", "3", "1", "2"))
)

colData <- data.frame(
  sample = colnames(filtered_counts),
  treatment = factor(c(rep("CK", 2), rep("MB", 2), rep("SCR", 2))), 
  replicate = factor(c("A", "B", "A", "B", "A", "B"))
)


# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = colData,
                              design = ~ treatment + replicate)

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)

# Convert results to data frame for plotting
res <- results(dds, contrast = c("treatment", "3wk_Ckm", "3wk_SCR"))
#res <- results(dds, contrast = c("treatment", "MB", "SCR"))
res_df <- as.data.frame(res)

# Check if res_df has the expected columns
if(!all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_df))) {
  stop("res_df does not contain expected columns")
}

# Define file path
file_path <- "C:/Users/harum/OneDrive/Documents/Rplot/HITI/Revision/3wkCKM_DES_1.pdf"

# Attempt to open the PDF device and ensure it is closed properly
pdf(file_path, width = 4, height = 4)
on.exit(dev.off(), add = TRUE)

# Create the volcano plot
with(res_df, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "3wkCkm vs SCR", xlim = c(-8.8, 8), col = "grey"))
with(subset(res_df, abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "lightblue"))
with(subset(res_df, padj < 0.05), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))

# Highlight specific genes
highlighted_genes <- c("Ckm", "Factor9")
highlighted_gene_index <- which(tolower(rownames(res_df)) %in% tolower(highlighted_genes))
if (length(highlighted_gene_index) > 0) {
  points(res_df$log2FoldChange[highlighted_gene_index], 
         -log10(res_df$pvalue)[highlighted_gene_index], 
         pch = 20, col = "red", cex = 1.5)
} else {
  cat("Highlighted genes are not found in the analysis.\n")
}

# Extract genes with adjusted p-value < 0.0005
significant_genes <- subset(res_df, padj < 0.0005)

# Annotate the significant genes on the volcano plot
if (nrow(significant_genes) > 0) {
  for (i in 1:nrow(significant_genes)) {
    text(x = significant_genes$log2FoldChange[i], y = -log10(significant_genes$pvalue[i]), 
         labels = rownames(significant_genes)[i], col = "black", adj = c(0, 0.5))
  }
} else {
  cat("No significant genes found.\n")
}

# Ensure the PDF device is closed
dev.off()

significant_genes_list <- rownames(res_df)[res_df$pvalue < 0.0005 & abs(res_df$log2FoldChange) > 1]
# Print significant genes
print(significant_genes_list)
