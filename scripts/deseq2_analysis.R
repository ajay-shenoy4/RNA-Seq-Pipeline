#!/usr/bin/env Rscript

# RNA-seq Differential Expression Analysis
# Supports both STAR+featureCounts and Salmon quantification

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

cat("=======================================================\n")
cat("RNA-seq Differential Expression Analysis\n")
cat("=======================================================\n\n")

# Create output directory
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. SALMON-BASED ANALYSIS (Primary Method)
# ============================================================================

cat("--- SALMON ANALYSIS ---\n\n")

# Find Salmon quant files
salmon_files <- list.files("results/salmon", 
                           pattern = "quant.sf$", 
                           recursive = TRUE, 
                           full.names = TRUE)

if (length(salmon_files) > 0) {
  cat(sprintf("Found %d Salmon quantification files\n", length(salmon_files)))
  
  # Create sample metadata
  sample_names <- gsub("_salmon/quant.sf", "", basename(dirname(salmon_files)))
  conditions <- ifelse(grepl("WT", sample_names), "WT", "KO")
  
  coldata_salmon <- data.frame(
    sample = sample_names,
    condition = factor(conditions, levels = c("WT", "KO")),
    files = salmon_files,
    row.names = sample_names
  )
  
  cat("\nSalmon sample metadata:\n")
  print(coldata_salmon[, c("sample", "condition")])
  
  # Create tx2gene mapping from GTF
  gtf_file <- "data/reference/gencode.v46.annotation.gtf"
  
  if (file.exists(gtf_file)) {
    cat("\nCreating transcript-to-gene mapping from GTF...\n")
    gtf <- read.table(gtf_file, sep = "\t", comment.char = "#", 
                      stringsAsFactors = FALSE)
    
    # Extract transcript and gene IDs
    tx2gene <- gtf[gtf$V3 == "transcript", ] %>%
      mutate(
        transcript_id = gsub('.*transcript_id "([^"]+)".*', '\\1', V9),
        gene_id = gsub('.*gene_id "([^"]+)".*', '\\1', V9)
      ) %>%
      select(transcript_id, gene_id) %>%
      distinct()
    
    cat(sprintf("Created tx2gene mapping with %d transcripts\n", nrow(tx2gene)))
    
    # Import Salmon counts with tximport
    cat("\nImporting Salmon counts...\n")
    txi <- tximport(coldata_salmon$files, 
                    type = "salmon", 
                    tx2gene = tx2gene,
                    ignoreTxVersion = TRUE)
    
    cat(sprintf("Imported counts for %d genes across %d samples\n", 
                nrow(txi$counts), ncol(txi$counts)))
    
    # Create DESeq2 object from Salmon
    dds_salmon <- DESeqDataSetFromTximport(
      txi,
      colData = coldata_salmon,
      design = ~ condition
    )
    
    # Filter low counts
    keep <- rowSums(counts(dds_salmon)) >= 10
    dds_salmon <- dds_salmon[keep, ]
    cat(sprintf("Kept %d genes after filtering (≥10 total reads)\n", sum(keep)))
    
    # Run DESeq2
    cat("\nRunning DESeq2 on Salmon data...\n")
    dds_salmon <- DESeq(dds_salmon)
    res_salmon <- results(dds_salmon, contrast = c("condition", "KO", "WT"))
    res_salmon <- res_salmon[order(res_salmon$pvalue), ]
    
    cat("\nSalmon DESeq2 summary:\n")
    summary(res_salmon)
    
    # Save Salmon results
    res_salmon_df <- as.data.frame(res_salmon)
    res_salmon_df$gene_id <- rownames(res_salmon_df)
    res_salmon_df <- res_salmon_df[, c("gene_id", "baseMean", "log2FoldChange", 
                                        "lfcSE", "stat", "pvalue", "padj")]
    write.csv(res_salmon_df, "results/deseq2/salmon_all_genes.csv", row.names = FALSE)
    
    sig_salmon <- subset(res_salmon_df, padj < 0.05)
    write.csv(sig_salmon, "results/deseq2/salmon_significant_genes.csv", row.names = FALSE)
    
    cat(sprintf("\nSalmon: Found %d significantly DE genes (padj < 0.05)\n", nrow(sig_salmon)))
    cat(sprintf("  - %d upregulated in KO\n", sum(sig_salmon$log2FoldChange > 0)))
    cat(sprintf("  - %d downregulated in KO\n", sum(sig_salmon$log2FoldChange < 0)))
    
    # Visualizations - Salmon
    cat("\nGenerating Salmon visualizations...\n")
    
    # MA Plot
    pdf("results/deseq2/salmon_ma_plot.pdf", width = 8, height = 6)
    plotMA(res_salmon, main = "Salmon: MA Plot (KO vs WT)", ylim = c(-5, 5), alpha = 0.05)
    dev.off()
    
    # Volcano Plot
    volcano_salmon <- as.data.frame(res_salmon) %>%
      mutate(
        significance = case_when(
          padj < 0.05 & log2FoldChange > 1 ~ "Up",
          padj < 0.05 & log2FoldChange < -1 ~ "Down",
          padj < 0.05 ~ "Sig",
          TRUE ~ "NS"
        )
      )
    
    p_volcano <- ggplot(volcano_salmon, 
                        aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(
        values = c("Up" = "red", "Down" = "blue", "Sig" = "orange", "NS" = "gray70")
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      labs(title = "Salmon: Volcano Plot (KO vs WT)",
           x = "log2(Fold Change)", y = "-log10(p-value)") +
      theme_bw()
    
    ggsave("results/deseq2/salmon_volcano_plot.pdf", p_volcano, width = 8, height = 6)
    
    # PCA Plot
    vsd_salmon <- vst(dds_salmon, blind = FALSE)
    
    pdf("results/deseq2/salmon_pca_plot.pdf", width = 7, height = 5)
    plotPCA(vsd_salmon, intgroup = "condition") + 
      theme_bw() +
      ggtitle("Salmon: PCA - Sample Clustering")
    dev.off()
    
    # Sample Distance Heatmap
    sampleDists <- dist(t(assay(vsd_salmon)))
    sampleDistMatrix <- as.matrix(sampleDists)
    colnames(sampleDistMatrix) <- NULL
    
    pdf("results/deseq2/salmon_sample_distance_heatmap.pdf", width = 7, height = 6)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
             main = "Salmon: Sample-to-Sample Distances")
    dev.off()
    
    # Top Genes Heatmap
    top_genes <- head(order(res_salmon$padj), 50)
    mat <- assay(vsd_salmon)[rownames(res_salmon)[top_genes], ]
    mat <- mat - rowMeans(mat)
    
    annotation_col <- data.frame(
      Condition = coldata_salmon$condition,
      row.names = rownames(coldata_salmon)
    )
    
    pdf("results/deseq2/salmon_top50_genes_heatmap.pdf", width = 8, height = 10)
    pheatmap(mat,
             annotation_col = annotation_col,
             cluster_rows = TRUE,
             show_rownames = FALSE,
             main = "Salmon: Top 50 DE Genes",
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
    
  } else {
    cat("GTF file not found, skipping tx2gene mapping\n")
  }
  
} else {
  cat("No Salmon quantification files found\n")
}

# ============================================================================
# 2. STAR+featureCounts ANALYSIS (If available)
# ============================================================================

cat("\n\n--- STAR+featureCounts ANALYSIS ---\n\n")

counts_file <- "results/counts/gene_counts.txt"

if (file.exists(counts_file)) {
  cat("Found featureCounts output\n")
  
  counts_data <- read.table(counts_file, header = TRUE, row.names = 1, 
                            skip = 1, check.names = FALSE)
  
  count_matrix <- as.matrix(counts_data[, 6:ncol(counts_data)])
  colnames(count_matrix) <- gsub(".*/(.*?)_Aligned.*", "\\1", colnames(count_matrix))
  
  sample_names <- colnames(count_matrix)
  conditions <- ifelse(grepl("WT", sample_names), "WT", "KO")
  
  coldata_star <- data.frame(
    sample = sample_names,
    condition = factor(conditions, levels = c("WT", "KO"))
  )
  rownames(coldata_star) <- sample_names
  
  cat("\nSTAR sample metadata:\n")
  print(coldata_star)
  
  dds_star <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata_star,
    design = ~ condition
  )
  
  keep <- rowSums(counts(dds_star)) >= 10
  dds_star <- dds_star[keep, ]
  
  cat("\nRunning DESeq2 on STAR data...\n")
  dds_star <- DESeq(dds_star)
  res_star <- results(dds_star, contrast = c("condition", "KO", "WT"))
  res_star <- res_star[order(res_star$pvalue), ]
  
  cat("\nSTAR DESeq2 summary:\n")
  summary(res_star)
  
  # Save STAR results
  res_star_df <- as.data.frame(res_star)
  res_star_df$gene_id <- rownames(res_star_df)
  res_star_df <- res_star_df[, c("gene_id", "baseMean", "log2FoldChange", 
                                  "lfcSE", "stat", "pvalue", "padj")]
  write.csv(res_star_df, "results/deseq2/star_all_genes.csv", row.names = FALSE)
  
  sig_star <- subset(res_star_df, padj < 0.05)
  write.csv(sig_star, "results/deseq2/star_significant_genes.csv", row.names = FALSE)
  
  cat(sprintf("\nSTAR: Found %d significantly DE genes (padj < 0.05)\n", nrow(sig_star)))
  
  # Comparison between methods if both exist
  if (exists("res_salmon") && exists("res_star")) {
    cat("\n\n--- COMPARING SALMON vs STAR ---\n\n")
    
    common_genes <- intersect(rownames(res_salmon), rownames(res_star))
    cat(sprintf("Common genes between methods: %d\n", length(common_genes)))
    
    comparison_df <- data.frame(
      gene_id = common_genes,
      salmon_log2FC = res_salmon[common_genes, "log2FoldChange"],
      star_log2FC = res_star[common_genes, "log2FoldChange"],
      salmon_padj = res_salmon[common_genes, "padj"],
      star_padj = res_star[common_genes, "padj"]
    )
    
    # Correlation plot
    p_compare <- ggplot(comparison_df, aes(x = salmon_log2FC, y = star_log2FC)) +
      geom_point(alpha = 0.3) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(title = "Salmon vs STAR: log2(Fold Change) Comparison",
           x = "Salmon log2FC", y = "STAR log2FC") +
      theme_bw()
    
    ggsave("results/deseq2/method_comparison.pdf", p_compare, width = 7, height = 7)
    
    cor_val <- cor(comparison_df$salmon_log2FC, comparison_df$star_log2FC, 
                   use = "complete.obs")
    cat(sprintf("Correlation between methods: %.3f\n", cor_val))
  }
  
} else {
  cat("No featureCounts output found (STAR workflow not run)\n")
}

# ============================================================================
# Save Session Info
# ============================================================================

writeLines(capture.output(sessionInfo()), 
           "results/deseq2/session_info.txt")

cat("\n=======================================================\n")
cat("✓ Analysis complete!\n")
cat("=======================================================\n")
cat("\nGenerated files in results/deseq2/:\n")
cat("  Salmon results:\n")
cat("    - salmon_all_genes.csv\n")
cat("    - salmon_significant_genes.csv\n")
cat("    - salmon_*.pdf (visualizations)\n")
if (file.exists(counts_file)) {
  cat("  STAR results:\n")
  cat("    - star_all_genes.csv\n")
  cat("    - star_significant_genes.csv\n")
  cat("    - method_comparison.pdf\n")
}
cat("\n")
