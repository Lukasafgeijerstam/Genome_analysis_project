# Load necessary libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Set working directory where count files are located
setwd("/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts_against_R7/count_files")

# List all count files 
count_files <- list.files(pattern = "*.txt")

# Read count files and extract counts
count_list <- lapply(count_files, function(file) {
  df <- read.delim(file, comment.char="#", header=TRUE, row.names=1)
  df <- df[, ncol(df), drop=FALSE]  
  colnames(df) <- gsub("\\.txt$", "", file)
  return(df)
})

# Combine into one count matrix
count_matrix <- do.call(cbind, count_list)

# Remove zero count genes
keep <- rowSums(count_matrix) > 0
count_matrix_filtered <- count_matrix[keep, ]

# Create sample metadata
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = c(rep("R7", 3), rep("HP126", 3))
)

# Convert condition to factor
sample_info$condition <- factor(sample_info$condition, levels = c("R7", "HP126"))

# Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Filter genes with less than 10 counts in at least 2 samples per subspecies
keep <- rowSums(counts(dds)[, dds$condition == "R7"]) > 0 & 
        rowSums(counts(dds)[, dds$condition == "HP126"]) > 0
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Get DE results (HP126 vs R7)
res <- results(dds, contrast = c("condition", "HP126", "R7"))

# Order results by adjusted p-value
res_ordered <- res[order(res$padj),]

# Save results to CSV
write.csv(as.data.frame(res_ordered), file = "DESeq2_results_HP126_vs_R7.csv")

# Plot MA plot
plotMA(res, ylim=c(-10,10))

# Plot PCA plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

# Plot heatmap of top 20 DEGs
top_genes <- rownames(res_ordered)[1:20]
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))

annotation_col <- data.frame(Condition = sample_info$condition)
rownames(annotation_col) <- rownames(sample_info)

pheatmap(mat_scaled, annotation_col = annotation_col,
         main = "Heatmap: Top 20 DEGs",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         fontsize_row = 8)