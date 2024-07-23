# Load necessary libraries
rm(list = ls())
cat("\014") 

library(DESeq2)
library(readxl)
library(ggplot2)
library(biomaRt)
library(GenomicRanges)
library(karyoploteR)

# Read raw counts data
raw_counts <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/rnaSEQ/rna_readcount_mouse_invivo.csv", header = TRUE, row.names = 1)

# Check for and handle NAs
if (any(is.na(raw_counts))) {
  raw_counts <- raw_counts[complete.cases(raw_counts), ]
}

# Print the first few rows
print(head(raw_counts))

# Filter out genes with low read counts
filtered_counts <- raw_counts[rowSums(raw_counts) > 10, ]

# Create colData with treatment information
colData <- data.frame(
  sample = colnames(filtered_counts),
  treatment = factor(c(rep("3wk_Ckm", 4), rep("3wk_SCR", 4), rep("8wk_Ckm", 4), rep("8wk_SCR", 3), rep("NT", 2))),
  #treatment = factor(c(rep("CK", 2), rep("SCR", 2))), 
  replicate = factor(c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "1", "2"))
  #replicate = factor(c("A", "B", "A", "B"))
)

# Just checking
print(colData)

# Set the levels
colData$treatment <- relevel(colData$treatment, ref = "8wk_SCR")

# Check the levels to confirm
print(levels(colData$treatment))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = colData,
                              design = ~ treatment + replicate)

# Check colData and design
print(colData)
print(dds@design)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Specify contrast
res <- results(dds, contrast = c("treatment", "8wk_Ckm", "8wk_SCR"))

# Convert results to a dataframe
res_df <- as.data.frame(res)

# Print the first few rows
print(head(res_df))

# Proceed with mapping
if (nrow(signif_genes_with_baseMean) > 0) {
  sigs <- signif_genes_with_baseMean
  
  # Annotate DEGs with chromosomal information using biomaRt
  ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  
  # Map gene symbols to Ensembl IDs
  gene_symbols <- rownames(sigs)
  
  if (length(gene_symbols) == 0) {
    stop("No gene symbols found in the results. Please check your DESeq2 analysis and input data.")
  }
  
  # Full Ensembl query
  gene_map <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
                    filters='external_gene_name', 
                    values=gene_symbols, 
                    mart=ensembl)
  
  # Merge DEGs with gene map to get Ensembl IDs
  deg_mapped <- merge(as.data.frame(sigs), gene_map, by.x="row.names", by.y="external_gene_name")
  
  # Check if we have matching Ensembl IDs
  if (nrow(deg_mapped) == 0) {
    stop("No matching Ensembl IDs were found. Please check your gene identifiers.")
  }
  
  # Annotate DEGs with chromosomal information using Ensembl IDs
  gene_annotations <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand'), 
                            filters='ensembl_gene_id', 
                            values=deg_mapped$ensembl_gene_id, 
                            mart=ensembl)
  
  # Merge DEGs with chromosomal information
  deg_annotated <- merge(deg_mapped, gene_annotations, by="ensembl_gene_id")
  
  # Create GRanges object for plotting
  gr <- GRanges(seqnames=Rle(deg_annotated$chromosome_name), 
                ranges=IRanges(start=deg_annotated$start_position, end=deg_annotated$end_position))
  
  new_seqnames <- paste0("chr", seqnames(gr))
  
  # Create a new GRanges object with updated chromosome names
  gr_updated <- GRanges(seqnames=Rle(new_seqnames),
                        ranges=IRanges(start=start(gr), end=end(gr)))
  
  # Assign colors based on regulation direction
  gr_colors <- rep("grey", length(gr_updated))  # Initialize with grey color
  gr_colors[deg_annotated$log2FoldChange > 0 & deg_annotated$padj < 0.05] <- "red"   # Upregulated genes in red
  gr_colors[deg_annotated$log2FoldChange < 0 & deg_annotated$padj < 0.05] <- "blue"  # Downregulated genes in blue
  
  # Plot chromosomal distribution using karyoploteR with updated GRanges and colors
  suppressWarnings({
    kp <- plotKaryotype(genome="mm10")
    kpPlotRegions(kp, data=gr_updated, col=gr_colors)
  })
  
} else {
  stop("No significant DEGs found with padj < 0.05 and baseMean > 50")
}
