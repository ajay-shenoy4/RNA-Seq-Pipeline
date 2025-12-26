## 1. Sample Clustering & Quality Control

PCA Plot (salmon_pca_plot.pdf)

What it is: A Principal Component Analysis that reduces thousands of gene expressions into two dimensions (PC1 and PC2).

What to look for: You want to see the three WT samples cluster together and the three KO samples cluster together.

Success Metric: If PC1 (the x-axis) separates your WT and KO groups, it means the "Condition" is the primary driver of difference in your data, rather than technical noise.

Sample Distance Heatmap (salmon_sample_distance_heatmap.pdf)

What it is: A "Euclidean distance" map showing how similar each sample is to every other sample.

What to look for: Darker blue blocks on the diagonal. You should see a clear 3x3 square for WT and a 3x3 square for KO. This confirms high reproducibility between your biological replicates.

## 2. Global Differential ExpressionVolcano Plot (salmon_volcano_plot.pdf)

What it is: A scatter plot comparing statistical significance ($-\log_{10} p$-value) against the magnitude of change ($\log_2$ Fold Change).

**How to read it:**
- Right Side (Red): Genes upregulated in the KO.
- Left Side (Blue): Genes downregulated in the KO.
- Top: The most statistically significant genes.

Your Result: With 5,155 DE genes, you should see a dense "V" shape. Look for DUS1L on the far left (downregulated).

MA Plot (salmon_ma_plot.pdf)

What it is: Plots log-intensity (Mean Expression) on the x-axis and log-fold change on the y-axis.

What to look for: This helps identify if fold-change estimates are biased by how lowly or highly expressed a gene is. Red dots represent genes with a padj < 0.05.

## 3. Top Gene Expression

Top 50 DE Genes Heatmap (salmon_top50_genes_heatmap.pdf)

What it is: A zoomed-in look at the 50 genes with the lowest $p$-values.

What to look for: A clear "binary" pattern where one block of genes is bright red in WT and blue in KO (and vice-versa).

The Goal: This provides a "snapshot" of the most reliable biomarkers of your Knockout condition.

## 4. Tabular Data Output

salmon_all_genes.csv: The master list of all 22,213 genes analyzed. Useful for custom searches.

salmon_significant_genes.csv: A filtered list containing only the 5,155 genes that passed the statistical threshold ($p_{adj} < 0.05$). This is the list you would typically use for downstream Gene Ontology (GO) or Pathway analysis.

## Summary of Biological Findings

**Based on the pipeline run, the DUS1L Knockout resulted in a massive transcriptomic shift:**

- Upregulated: ~2,533 genes
- Downregulated: ~2,622 genes
- Validation: The pipeline correctly identified the targeted depletion of the DUS1L transcript, confirming the experimental validity of the dataset.






















