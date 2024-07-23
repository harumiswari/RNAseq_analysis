# Load necessary libraries
library(DESeq2)
library(readxl)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Read raw counts data
raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/Raw_Read_Mouse.csv", header = TRUE, row.names = 1)
#raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/rna_readcount_mouse_invivo.csv", header = TRUE, row.names = 1)

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

#colData <- data.frame(
#  sample = colnames(filtered_counts),
#  treatment = factor(c(rep("3wk_Ckm", 4), rep("3wk_SCR", 4), rep("8wk_Ckm", 4), rep("8wk_SCR", 3), rep("NT", 2))),
#  replicate = factor(c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "1", "2"))
#)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = colData, design = ~ treatment + replicate)

# Run DESeq2 analysis
dds <- DESeq(dds)

#Perform differential expression analysis for  vs SCR
res <- results(dds, contrast = c("treatment", "CK", "SCR"))
res_df <- as.data.frame(res)

# Filter significant genes for GO analysis
sigs <- res_df[res_df$padj < 0.05 & res_df$baseMean > 50, ]

# Extract gene symbols for GO analysis
genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5, ])

# Perform GO enrichment analysis
GO_results <- enrichGO(gene = genes_to_test, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")


# Plot the top 10 GO terms based on the number of genes involved
library(ggplot2)
# Ensure 'GO_results' is converted to a data frame
GO_results_df <- as.data.frame(GO_results)

# Print all GO enrichment results
print(GO_results_df)

# Print only selected columns, for example: ID, Description, p.adjust, and gene count
print(GO_results_df[, c("ID", "Description", "p.adjust", "Count")])

# Sort results by adjusted p-value to see the most significant terms first
GO_results_sorted <- GO_results_df[order(GO_results_df$p.adjust), ]
print(GO_results_sorted)


top_GO_terms <- head(GO_results_sorted, 10)  # Adjust number as necessary
#top_GO_terms <- GO_results_sorted

# Create the plot
p <- ggplot(top_GO_terms, aes(x = reorder(Description, -Count), y = Count, fill = -log10(p.adjust))) +
  geom_col() +
  coord_flip() +  # Flips the x and y axes
  scale_fill_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
  labs(title = "Top 10 specific GO Terms by Gene Count", x = "GO Term", y = "Gene Count") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # Adjust margin for more spacing
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")  # Increase plot margin
  )

# Save the plot as a PDF file
ggsave("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/CK_BP_Top_GO_Terms.pdf", plot = p, width = 12, height = 8, units = "in", dpi = 300)
# Extend the list of keywords to cover more cancer-related processes
#keywords <- "cancer|tumorigenesis|oncogenesis|cell proliferation|apoptosis|angiogenesis|cell migration|DNA damage response|epithelial to mesenchymal transition"
#keywords <- "glycolysis|immune response|evasion|inflammatory response|DNA repair"
#keywords <- "positive regulation of cell division|negative regulation of cell cycle arrest|negative regulation of signal transduction by p53 class mediator|negative regulation of apoptotic process|anti-apoptosis|telomere maintenance|regulation of telomerase activity|positive regulation of angiogenesis|cell invasion|activation of immune response|response to DNA damage stimulus| chromosomal instability|immune evasion|epithelial-to-mesenchymal transition"
# Search for GO terms using the extended list of keywords
keywords <- "muscle contraction|muscle tissue regeneration"
cancer_related_go_terms <- grep(keywords, GO_results_df$Description, ignore.case = TRUE, value = TRUE)

# Print the matching GO terms
print(GO_results_df[GO_results_df$Description %in% cancer_related_go_terms, ])

cancer_GO_data <- GO_results_df[GO_results_df$Description %in% cancer_related_go_terms, ]

write.csv(cancer_GO_data, "C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/MB_BP2_Cancer_Related_GO_Terms.csv", row.names = FALSE)
# Convert GO_results to a data frame 
GO_results_df <- as.data.frame(GO_results)
GO_results_df$p.adjust <- as.numeric(GO_results_df$p.adjust)


# Check if GO_results_df contains any results
if (nrow(GO_results_df) > 0) {
  # Create the barplot using ggplot2 with gradient fill
  p <- ggplot(GO_results_df, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "blue", high = "red", name = "Adjusted p-value") +
    labs(title = "GO Biological Process Enrichment",
         x = "GO Term",
         y = "Gene Count") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),  # Adjust margin for more spacing
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold", margin = margin(b = 10))) # Margin for plot title
  
  # Print the plot
  print(p)
} else {
  print("No significant GO terms found.")
}


############## for heatmap################
read_data <- read.csv ("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/MBvsSCR_GO.csv")
# Load necessary library
library(pheatmap)

# Suppose your CSV file is named 'GO_data.csv' and located in the current working directory
GO_data <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/MBvsSCR_GO.csv")

# Prepare a matrix for the heatmap. Here, we focus on 'p.adjust' and 'Count'
heatmap_data <- as.matrix(GO_data[, c("p.adjust", "Count")])
rownames(heatmap_data) <- GO_data$Description  # Set GO term descriptions as row names

# Normalize p.adjust by converting to -log10 scale for better visualization
heatmap_data[, "p.adjust"] <- -log10(heatmap_data[, "p.adjust"])

# Generate the heatmap
pheatmap(heatmap_data,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row",  # Normalize data per row
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color scale
         main = "GO Term Enrichment Heatmap",
         legend_title = "Metrics")
