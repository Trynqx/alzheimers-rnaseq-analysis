# Alzheimer's RNA-seq analysis pipeline
# Dataset: GSE53697
# Author: Surya Kumar

# =========================
# 1. Load packages
# =========================
library(DESeq2)
library(tidyverse)
library(pheatmap)

# =========================
# 2. Create output folders if missing
# =========================
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

# =========================
# 3. Load data
# =========================
counts <- read.delim(
  "data/GSE53697_RNAseq_AD.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# =========================
# 4. Clean gene names
# =========================
counts <- counts[counts$GeneSymbol != "" & !is.na(counts$GeneSymbol), ]
counts$GeneSymbol <- make.unique(as.character(counts$GeneSymbol))

rownames(counts) <- counts$GeneSymbol
counts <- counts[, -1]

# =========================
# 5. Keep only raw counts
# =========================
count_matrix <- as.matrix(counts)
storage.mode(count_matrix) <- "integer"
count_matrix <- count_matrix[, grepl("_raw", colnames(count_matrix))]

# =========================
# 6. Create metadata
# =========================
sample_names <- colnames(count_matrix)

metadata <- data.frame(
  sample = sample_names,
  condition = ifelse(grepl("^C", sample_names), "Control", "AD")
)

rownames(metadata) <- metadata$sample
metadata$condition <- factor(metadata$condition, levels = c("Control", "AD"))

# Save metadata
write.csv(metadata, "results/metadata.csv", row.names = TRUE)

# =========================
# 7. Run DESeq2
# =========================
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

# =========================
# 8. Extract results
# =========================
res <- results(dds, contrast = c("condition", "AD", "Control"))
res <- as.data.frame(res)
res$gene <- rownames(res)

# Order by p-value
res <- res[order(res$pvalue), ]

# Save all results
write.csv(res, "results/all_deseq2_results.csv", row.names = FALSE)

# =========================
# 9. Exploratory significance
# =========================
res_p <- subset(res, !is.na(pvalue) & pvalue < 0.05)

sig_res <- subset(
  res,
  !is.na(pvalue) &
    pvalue < 0.05 &
    abs(log2FoldChange) > 1
)

write.csv(res_p, "results/pvalue_lt_0.05_results.csv", row.names = FALSE)
write.csv(sig_res, "results/exploratory_significant_genes.csv", row.names = FALSE)

cat("Genes with p < 0.05:", nrow(res_p), "\n")
cat("Genes with p < 0.05 and |log2FC| > 1:", nrow(sig_res), "\n")

# =========================
# 10. Volcano plot
# =========================
res$significant <- ifelse(
  !is.na(res$pvalue) & res$pvalue < 0.05 & abs(res$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

res_plot <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]

volcano_plot <- ggplot(
  res_plot,
  aes(x = log2FoldChange, y = -log10(pvalue), color = significant)
) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  ggtitle("Volcano Plot: Alzheimer's vs Control") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red"))

print(volcano_plot)
ggsave("figures/volcano.png", plot = volcano_plot, width = 7, height = 5)

# =========================
# 11. PCA plot
# =========================
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("PCA: Alzheimer's vs Control") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

print(pca_plot)
ggsave("figures/PCA.png", plot = pca_plot, width = 7, height = 5)

# =========================
# 12. Heatmap of top 30 genes
# =========================
top_n <- min(30, nrow(sig_res))

if (top_n > 1) {
  top_genes <- sig_res$gene[1:top_n]
  
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  
  annotation_df <- data.frame(condition = metadata$condition)
  rownames(annotation_df) <- rownames(metadata)
  
  pheatmap(
    mat,
    annotation_col = annotation_df,
    filename = "figures/heatmap.png"
  )
}

# =========================
# 13. Save normalized counts
# =========================
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts), "results/normalized_counts.csv")

# =========================
# 14. Session summary
# =========================
cat("\nAnalysis complete.\n")
cat("Output files saved in 'results/' and 'figures/' folders.\n")