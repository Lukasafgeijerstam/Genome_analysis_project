# DE_analysis.R

# Load libraries
library(DESeq2)

# --- Step 1: Read count data ---
# Adjust paths to your files
r7_counts_file_62 <- "/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts/R7_counts_62.txt"
r7_counts_file_63 <- "/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts/R7_counts_63.txt"
r7_counts_file_64 <- "/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts/R7_counts_64.txt"
hp126_counts_file_59 <- "/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts/HP126_counts_59.txt"
hp126_counts_file_60 <- "/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts/HP126_counts_60.txt"
hp126_counts_file_61 <- "/home/lukasa/Genome_analysis_project/08-read_counting/featureCounts/HP126_counts_61.txt"

# Read counts as tables, skipping featureCounts comments 
r7_62 <- read.table(r7_counts_file_62, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)
r7_63 <- read.table(r7_counts_file_63, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)
r7_64 <- read.table(r7_counts_file_64, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)
hp126_59 <- read.table(hp126_counts_file_59, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)
hp126_60 <- read.table(hp126_counts_file_60, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)
hp126_61 <- read.table(hp126_counts_file_61, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)

# Extract the counts column 
r7_62_counts <- r7_62[, ncol(r7_62), drop = FALSE]
r7_63_counts <- r7_63[, ncol(r7_63), drop = FALSE]
r7_64_counts <- r7_64[, ncol(r7_64), drop = FALSE]
hp126_59_counts <- hp126_59[, ncol(hp126_59), drop = FALSE]
hp126_60_counts <- hp126_60[, ncol(hp126_60), drop = FALSE]
hp126_61_counts <- hp126_61[, ncol(hp126_61), drop = FALSE]

# Combine into one count matrix
count_matrix <- cbind(R7_62 = r7_62_counts[,1], R7_63 = r7_63_counts[,1], R7_64 = r7_64_counts[,1], 
HP126_59 = hp126_59_counts[,1], HP126_60 = hp126_60_counts[,1], HP126_61 = hp126_61_counts[,1])

# --- Step 2: Create sample metadata ---
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c("R7", "R7", "R7", "HP", "HP126", "HP126")),
)

# --- Step 3: Create DESeqDataSet ---
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Pre-filtering: remove genes with zero counts across all samples
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# --- Step 4: Run DESeq ---
dds <- DESeq(dds)

# --- Step 5: Get results ---
res <- results(dds)
res <- res[order(res$padj), ]  # Sort by adjusted p-value

# --- Step 6: Save results ---
write.csv(as.data.frame(res), file = "DE_results_R7_vs_HP126.csv")

# --- Step 7: Basic plots ---

# MA-plot
png("MAplot_DESeq2.png")
plotMA(res, main="DESeq2 MA-plot", ylim=c(-5,5))
dev.off()

# Sample distance heatmap and PCA require replicates, so omitted here

# Print summary to console
print(summary(res))



