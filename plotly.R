# Load necessary libraries
library(DESeq2)
library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(plotly)

# Read raw counts data
raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/Raw_Read_Mouse.csv", header = TRUE, row.names = 1)
# raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/rna_readcount_mouse_invivo.csv", header = TRUE, row.names = 1)

# Check for and handle NAs
if (any(is.na(raw_counts))) {
  raw_counts <- na.omit(raw_counts)
}

# Filter out genes with low read counts
filtered_counts <- raw_counts[rowSums(raw_counts) > 10, ]

# Create colData with treatment information
colData <- data.frame(
  sample = colnames(filtered_counts),
  treatment = factor(c(rep("CK", 2), rep("MB", 2), rep("SCR", 2))),
  replicate = factor(rep(c("A", "B"), 3))
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = colData, design = ~ treatment + replicate)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Perform differential expression analysis for MB vs SCR
res <- results(dds, contrast = c("treatment", "MB", "SCR"))
res_df <- as.data.frame(res)

# Filter significant genes for GO analysis
sigs <- res_df[res_df$padj < 0.05 & res_df$baseMean > 50, ]

# Extract gene symbols for GO analysis
genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5, ])

# Perform GO enrichment analysis for BP, MF, and CC
GO_results_BP <- enrichGO(gene = genes_to_test, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
GO_results_MF <- enrichGO(gene = genes_to_test, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "MF")
GO_results_CC <- enrichGO(gene = genes_to_test, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC")

# Convert results to data frames
GO_results_BP_df <- as.data.frame(GO_results_BP)
GO_results_MF_df <- as.data.frame(GO_results_MF)
GO_results_CC_df <- as.data.frame(GO_results_CC)

# Add a column to indicate the ontology
GO_results_BP_df$Category <- "BP"
GO_results_MF_df$Category <- "MF"
GO_results_CC_df$Category <- "CC"

# Combine the results into one data frame
GO_results_combined <- bind_rows(GO_results_BP_df, GO_results_MF_df, GO_results_CC_df)

# Sort results by adjusted p-value
GO_results_sorted <- GO_results_combined %>%
  arrange(p.adjust) %>%
  group_by(Category) %>%
  top_n(-30, p.adjust) %>%
  ungroup()

# Plot with Plotly
library(plotly)

fig <- plot_ly(GO_results_sorted, x = ~-log10(p.adjust), y = ~Description, color = ~Category, size = ~Count,
               type = 'scatter', mode = 'markers', text = ~paste("GO ID:", ID, "<br>Adjusted P-value:", p.adjust, "<br>Gene Count:", Count))

fig <- fig %>% layout(title = 'Top Enriched GO Terms',
                      xaxis = list(title = '-log10(Adjusted P-value)'),
                      yaxis = list(title = 'GO Term'),
                      legend = list(title = list(text = 'Category')))

fig
