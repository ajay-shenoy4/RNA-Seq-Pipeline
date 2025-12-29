#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggrepel)
})

cat("=======================================================\n")
cat("RNA-seq Differential Expression Analysis\n")
cat("=======================================================\n\n")

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

salmon_files <- list.files("results/salmon", pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)

if (length(salmon_files) > 0) {
  sample_names <- gsub("_salmon/quant.sf", "", basename(dirname(salmon_files)))
  conditions <- ifelse(grepl("WT", sample_names), "WT", "KO")
  coldata_salmon <- data.frame(sample = sample_names, condition = factor(conditions, levels = c("WT", "KO")), files = salmon_files, row.names = sample_names)

  gtf_file <- "data/reference/gencode.v46.annotation.gtf"
  if (file.exists(gtf_file)) {
    gtf_lines <- readLines(gtf_file)
    gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)] 
    tx2gene <- do.call(rbind, Filter(Negate(is.null), lapply(gtf_lines, function(line) {
      if (grepl("\ttranscript\t", line)) {
        tx_id <- gsub('.*transcript_id "([^"]+)".*', "\\1", line)
        gene_id <- gsub('.*gene_id "([^"]+)".*', "\\1", line)
        # We will keep both the Ensembl ID and the Symbol if possible
        gene_name <- gsub('.*gene_name "([^"]+)".*', "\\1", line)
        return(c(sub("\\..*", "", tx_id), sub("\\..*", "", gene_id), gene_name))
      }
      return(NULL)
    })))
    tx2gene <- as.data.frame(tx2gene, stringsAsFactors = FALSE)
    colnames(tx2gene) <- c("transcript_id", "gene_id", "gene_name")
    
    txi <- tximport(coldata_salmon$files, type = "salmon", tx2gene = tx2gene[,1:2], ignoreTxVersion = TRUE)
    dds <- DESeqDataSetFromTximport(txi, colData = coldata_salmon, design = ~condition)
    dds <- DESeq(dds[rowSums(counts(dds)) >= 10, ])
    res <- results(dds, contrast = c("condition", "KO", "WT"))
    res <- res[order(res$pvalue), ]

    # Map Ensembl IDs to Symbols for Plotting
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df$symbol <- tx2gene$gene_name[match(res_df$gene_id, tx2gene$gene_id)]
    
    write.csv(res_df, "results/deseq2/salmon_all_genes.csv", row.names = FALSE)

    # 1. Volcano Plot with proper labels
    volcano_data <- res_df %>%
      mutate(neglog10p = -log10(pvalue),
             sig = case_when(padj < 0.05 & log2FoldChange > 1 ~ "Up",
                             padj < 0.05 & log2FoldChange < -1 ~ "Down",
                             padj < 0.05 ~ "Sig", TRUE ~ "NS"),
             label = ifelse(neglog10p > 20 | symbol == "DUS1L", symbol, NA))

    p_volc <- ggplot(volcano_data, aes(x=log2FoldChange, y=neglog10p, color=sig)) +
      geom_point(alpha=0.4) +
      geom_text_repel(aes(label=label), max.overlaps = 15) +
      scale_color_manual(values=c("Up"="red", "Down"="blue", "Sig"="orange", "NS"="grey")) +
      theme_bw() + labs(title="Volcano Plot: KO vs WT")
    ggsave("results/deseq2/salmon_volcano_plot.pdf", p_volc, width=8, height=6)

    # 2. Functional Heatmap with Error Handling
    vsd <- vst(dds, blind=FALSE)
    top_idx <- head(order(res$padj), 50)
    mat <- assay(vsd)[top_idx, ]
    rownames(mat) <- res_df$symbol[match(rownames(mat), res_df$gene_id)]
    mat <- mat - rowMeans(mat)

    # Attempt GO Enrichment
    ego <- tryCatch({
      enrichGO(gene = rownames(mat), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
    }, error = function(e) NULL)

    ann_row <- NULL
    if (!is.null(ego) && nrow(ego) > 0) {
      gene2go <- setNames(ego@result$Description[match(rownames(mat), ego@result$geneID)], rownames(mat))
      gene2go[is.na(gene2go)] <- "Other"
      ann_row <- data.frame(Process = gene2go, row.names = rownames(mat))
    }

    pdf("results/deseq2/salmon_top50_genes_heatmap.pdf", width = 10, height = 12)
    pheatmap(mat, annotation_col = data.frame(Cond = coldata_salmon$condition, row.names = colnames(mat)),
             annotation_row = ann_row, main = "Top 50 DE Genes", color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
    
    # 3. Standard Plots
    pdf("results/deseq2/salmon_pca_plot.pdf", width=7, height=5)
    print(plotPCA(vsd, intgroup="condition") + theme_bw())
    dev.off()
  }
}
cat("\n✓ Analysis complete!\n")
