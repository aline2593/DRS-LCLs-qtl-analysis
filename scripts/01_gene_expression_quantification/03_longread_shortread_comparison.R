#!/usr/bin/env Rscript
################################################################################
# Long-read vs Short-read Gene Expression Comparison
################################################################################
# Description: Compare gene expression quantified from Direct RNA-seq (nanopore)
#              versus standard short-read RNA-seq (Illumina) in the same samples
              (in this study are 60 LCLs from GEUVADIS from European individuals).
#              Analyzes correlation, overlap, and characteristics of detected genes.
#
# Input:  - Filtered nanopore expression BED (from script 02)
#         - Filtered Illumina expression BED (from Delaneau et al. 2019, Science;
            50% expression re-calculated on the subset of 60 LCLs common between
            the two studies)
#         - Gene length annotation
#
# Output: - Correlation plots and statistics
#         - Expression level comparisons
#         - Gene length comparisons
#         - Overlap analysis
#
# Dependencies: R packages - data.table, Hmisc, corrplot, ggplot2, ggsignif, tidyr
################################################################################

library(data.table)
library(Hmisc)
library(corrplot)
library(ggplot2)
library(ggsignif)
library(tidyr)

# Set working directory
setwd("/path/to/expression/comparison/")

################################################################################
# STEP 1: LOAD DATA
################################################################################

cat("Loading expression data...\n")

# Load filtered expression data (50% expression threshold applied)
# Both files should have same structure: #chr | start | end | gene_id | . | strand | samples...
nanopore_bed <- fread("gene_expression_filtered50_nanopore.bed.gz", data.table = FALSE)
illumina_bed <- fread("gene_expression_filtered50_illumina.bed.gz", data.table = FALSE)

cat("Nanopore genes:", nrow(nanopore_bed), "\n")
cat("Illumina genes:", nrow(illumina_bed), "\n")

################################################################################
# STEP 2: STANDARDIZE GENE IDs
################################################################################
# Remove version numbers from Ensembl IDs (e.g., ENSG00000223972.5 -> ENSG00000223972)

nanopore_bed$id <- sub("\\..*", "", nanopore_bed$id)
illumina_bed$gene <- sub("\\..*", "", illumina_bed$gene)

################################################################################
# STEP 3: IDENTIFY OVERLAPPING GENES
################################################################################
# Find genes quantified in both technologies

overlapping_genes <- intersect(nanopore_bed$id, illumina_bed$gene)
cat("\nGenes detected in both technologies:", length(overlapping_genes), "\n")

# Extract overlapping genes and remove BED coordinate columns
nanopore_overlap <- nanopore_bed[nanopore_bed$id %in% overlapping_genes, -c(1,2,3,5,6)]
illumina_overlap <- illumina_bed[illumina_bed$gene %in% overlapping_genes, -c(1,2,3,5,6)]

# Match gene order between datasets
nanopore_overlap <- nanopore_overlap[match(illumina_overlap$gene, nanopore_overlap$id), ]

# Verify gene order matches
stopifnot(sum(nanopore_overlap$id != illumina_overlap$gene) == 0)
cat("Gene order verified: PASS\n")

# Match sample order (same samples should be in same order)
illumina_overlap <- illumina_overlap[, match(colnames(nanopore_overlap), colnames(illumina_overlap))]
stopifnot(sum(colnames(nanopore_overlap) != colnames(illumina_overlap)) == 0)
cat("Sample order verified: PASS\n\n")

################################################################################
# STEP 4: CALCULATE CORRELATION BETWEEN TECHNOLOGIES
################################################################################
# Compare gene expression correlation across samples and within technologies

cat("Calculating correlations...\n")

# Prefix column names to distinguish technologies
colnames(nanopore_overlap)[-1] <- paste0("N_", colnames(nanopore_overlap)[-1])
colnames(illumina_overlap)[-1] <- paste0("I_", colnames(illumina_overlap)[-1])

# Combine datasets for correlation analysis
combined_data <- cbind(nanopore_overlap[, -1], illumina_overlap[, -1])

# Calculate Spearman correlation (robust to non-linear relationships)
# Structure: first 60 cols = nanopore, next 60 cols = Illumina
correlation_matrix <- rcorr(as.matrix(combined_data), type = "spearman")

# Save full correlation heatmap
pdf("correlation_heatmap_nanopore_vs_illumina.pdf", width = 12, height = 12)
corrplot(correlation_matrix$r,
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.6)
dev.off()

################################################################################
# STEP 5: ANALYZE MATCHED SAMPLE CORRELATIONS
################################################################################
# For each sample, correlation between nanopore and Illumina quantification

# Extract cross-technology correlations (nanopore rows vs Illumina columns)
cross_tech_corr <- correlation_matrix$r[1:60, 61:120]

# Get diagonal = matched sample correlations (sample1_nano vs sample1_illumina)
matched_correlations <- diag(cross_tech_corr)

cat("Matched sample correlation statistics:\n")
cat("  Median:", median(matched_correlations), "\n")
cat("  Mean:", mean(matched_correlations), "\n")
cat("  Range:", range(matched_correlations), "\n\n")

# Plot histogram of matched correlations
pdf("matched_sample_correlations.pdf", width = 8, height = 6)
hist(matched_correlations,
     main = "Gene Expression Correlation\nNanopore vs Illumina (Matched Samples)",
     xlab = "Spearman correlation coefficient",
     col = "cadetblue4",
     border = "white",
     breaks = 20,
     xlim = c(0, 1))
abline(v = median(matched_correlations), col = "red", lwd = 3, lty = 2)
text(x = median(matched_correlations),
     y = par("usr")[4] * 0.9,
     labels = paste("Median =", round(median(matched_correlations), 3)),
     pos = 4, col = "red")
dev.off()

################################################################################
# STEP 6: ANALYZE ALL PAIRWISE CORRELATIONS
################################################################################
# Include all sample pairs (not just matched)

all_pairs <- c(cross_tech_corr[lower.tri(cross_tech_corr)], matched_correlations)

pdf("all_pairs_correlations.pdf", width = 8, height = 6)
hist(all_pairs,
     main = "All Pairwise Correlations\nNanopore vs Illumina",
     xlab = "Spearman correlation coefficient",
     col = "grey",
     border = "white",
     breaks = 20,
     xlim = c(0, 1))
abline(v = median(cross_tech_corr), col = "red", lwd = 3, lty = 2)
text(x = median(cross_tech_corr),
     y = par("usr")[4] * 0.9,
     labels = paste("Median =", round(median(cross_tech_corr), 3)),
     pos = 4, col = "red")
dev.off()

################################################################################
# STEP 7: WITHIN-TECHNOLOGY CORRELATIONS
################################################################################

# Illumina vs Illumina
illumina_self_corr <- correlation_matrix$r[61:120, 61:120]
cat("Illumina self-correlation median:", median(illumina_self_corr[lower.tri(illumina_self_corr)]), "\n")

pdf("illumina_within_technology_correlation.pdf", width = 8, height = 6)
hist(illumina_self_corr[lower.tri(illumina_self_corr)],
     main = "Within-Technology Correlation\nIllumina",
     xlab = "Spearman correlation coefficient",
     col = "grey",
     border = "white",
     breaks = 20,
     xlim = c(0, 1))
abline(v = median(illumina_self_corr[lower.tri(illumina_self_corr)]),
       col = "red", lwd = 3, lty = 2)
dev.off()

# Nanopore vs Nanopore
nanopore_self_corr <- correlation_matrix$r[1:60, 1:60]
cat("Nanopore self-correlation median:", median(nanopore_self_corr[lower.tri(nanopore_self_corr)]), "\n\n")

pdf("nanopore_within_technology_correlation.pdf", width = 8, height = 6)
hist(nanopore_self_corr[lower.tri(nanopore_self_corr)],
     main = "Within-Technology Correlation\nDirect RNA-seq",
     xlab = "Spearman correlation coefficient",
     col = "grey",
     border = "white",
     breaks = 20,
     xlim = c(0, 1))
abline(v = median(nanopore_self_corr[lower.tri(nanopore_self_corr)]),
       col = "red", lwd = 3, lty = 2)
dev.off()

################################################################################
# STEP 8: GENE-LEVEL EXPRESSION COMPARISON
################################################################################
# Compare expression levels of genes detected in both vs only one technology

cat("Analyzing gene expression characteristics...\n")

# Identify technology-specific genes
nanopore_only_genes <- nanopore_bed[!(nanopore_bed$id %in% overlapping_genes), ]
illumina_only_genes <- illumina_bed[!(illumina_bed$gene %in% overlapping_genes), ]

cat("Nanopore-only genes:", nrow(nanopore_only_genes), "\n")
cat("Illumina-only genes:", nrow(illumina_only_genes), "\n\n")

# Calculate mean expression per gene (across samples)
calc_mean_expression <- function(bed_data) {
  expr_cols <- bed_data[, 7:ncol(bed_data)]
  data.frame(mean_expr = rowMeans(expr_cols))
}

# Create datasets with mean expression and group labels
nanopore_all_mean <- calc_mean_expression(nanopore_bed)
nanopore_all_mean$Group <- "dRNA-seq (all)"

nanopore_only_mean <- calc_mean_expression(nanopore_only_genes)
nanopore_only_mean$Group <- "dRNA-seq only"

illumina_all_mean <- calc_mean_expression(illumina_bed)
illumina_all_mean$Group <- "Short-reads (all)"

illumina_only_mean <- calc_mean_expression(illumina_only_genes)
illumina_only_mean$Group <- "Short-reads only"

# Statistical tests
pval_illumina <- wilcox.test(illumina_all_mean$mean_expr,
                             illumina_only_mean$mean_expr,
                             alternative = "two.sided")$p.value

pval_nanopore <- wilcox.test(nanopore_all_mean$mean_expr,
                             nanopore_only_mean$mean_expr,
                             alternative = "two.sided")$p.value

cat("Statistical tests (Wilcoxon rank-sum):\n")
cat("  Illumina all vs only:", format(pval_illumina, scientific = TRUE), "\n")
cat("  Nanopore all vs only:", format(pval_nanopore, scientific = TRUE), "\n\n")

# Combine for plotting
combined_expr <- rbind(
  data.frame(mean = nanopore_only_mean$mean_expr, Group = "dRNA-seq only"),
  data.frame(mean = illumina_all_mean$mean_expr, Group = "Short-reads"),
  data.frame(mean = nanopore_all_mean$mean_expr, Group = "dRNA-seq"),
  data.frame(mean = illumina_only_mean$mean_expr, Group = "Short-reads only")
)

# Order factor levels for plotting
combined_expr$Group <- factor(combined_expr$Group,
                              levels = c("dRNA-seq", "dRNA-seq only",
                                       "Short-reads", "Short-reads only"))

# Create violin plot

expr_plot <- ggplot(combined_expr, aes(x = Group, y = mean, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_y_continuous(trans = 'log10') +
  theme_light(base_size = 14) +
  xlab("Gene Groups") +
  ylab("Mean Gene Expression (log10 RPKM)") +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold", vjust = 1.8),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("gene_expression_comparison_violin.pdf",
       plot = expr_plot,
       height = 8,
       width = 12)

################################################################################
# STEP 9: GENE LENGTH COMPARISON
################################################################################
# Compare gene lengths between technology-specific and shared genes

cat("Analyzing gene length characteristics...\n")

# Load gene length annotations
gene_lengths <- fread("/path/to/annotation/Gene_length.txt", data.table = FALSE)

# For Illumina data with GTF attributes format
illumina_lengths <- fread("/path/to/annotation/illumina_gene_lengths.txt",
                         data.table = FALSE)

# Parse length from GTF attributes if needed (format: gene_id=X;length=Y)
if ("gid" %in% colnames(illumina_lengths)) {
  illumina_lengths <- illumina_lengths %>%
    separate(gid, c("L", "length"), "=") %>%
    select(gene, length)
  illumina_lengths$length <- as.numeric(illumina_lengths$length)
}

# Standardize gene IDs
gene_lengths$Geneid <- sub("\\..*", "", gene_lengths$Geneid)
if ("gene" %in% colnames(illumina_lengths)) {
  illumina_lengths$gene <- sub("\\..*", "", illumina_lengths$gene)
}

# Merge expression data with gene lengths
nanopore_with_length <- merge(nanopore_bed, gene_lengths,
                              by.x = "id", by.y = "Geneid")
nanopore_only_with_length <- merge(nanopore_only_genes, gene_lengths,
                                   by.x = "id", by.y = "Geneid")

illumina_with_length <- merge(illumina_bed, illumina_lengths, by = "gene")
illumina_only_with_length <- merge(illumina_only_genes, illumina_lengths,
                                   by = "gene")

# Prepare datasets for plotting
prepare_length_data <- function(data, group_name, length_col) {
  data.frame(
    length = data[[length_col]],
    Group = group_name
  )
}

length_nanopore <- prepare_length_data(nanopore_with_length, "dRNA-seq", "Length")
length_nanopore_only <- prepare_length_data(nanopore_only_with_length,
                                           "dRNA-seq only", "Length")
length_illumina <- prepare_length_data(illumina_with_length, "Short-reads", "length")
length_illumina_only <- prepare_length_data(illumina_only_with_length,
                                           "Short-reads only", "length")

# Statistical tests
pval_length_illumina <- wilcox.test(length_illumina$length,
                                   length_illumina_only$length,
                                   alternative = "two.sided")$p.value

pval_length_nanopore <- wilcox.test(length_nanopore$length,
                                   length_nanopore_only$length,
                                   alternative = "two.sided")$p.value

cat("Gene length comparison (Wilcoxon rank-sum):\n")
cat("  Illumina all vs only:", format(pval_length_illumina, scientific = TRUE), "\n")
cat("  Nanopore all vs only:", format(pval_length_nanopore, scientific = TRUE), "\n\n")

# Combine for plotting
combined_length <- rbind(length_illumina, length_nanopore,
                        length_illumina_only, length_nanopore_only)

combined_length$Group <- factor(combined_length$Group,
                               levels = c("dRNA-seq", "dRNA-seq only",
                                        "Short-reads", "Short-reads only"))

# Create violin plot

length_plot <- ggplot(combined_length, aes(x = Group, y = length, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_y_continuous(trans = 'log10') +
  theme_light(base_size = 14) +
  xlab("Gene Groups") +
  ylab("Gene Length (bp, log10)") +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold", vjust = 1.8),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("gene_length_comparison_violin.pdf",
       plot = length_plot,
       height = 8,
       width = 12)

################################################################################
# STEP 10: SUMMARY STATISTICS
################################################################################

cat("\n=== SUMMARY STATISTICS ===\n\n")

cat("Gene Detection:\n")
cat("  Total genes (nanopore):", nrow(nanopore_bed), "\n")
cat("  Total genes (Illumina):", nrow(illumina_bed), "\n")
cat("  Overlapping genes:", length(overlapping_genes), "\n")
cat("  Nanopore-specific:", nrow(nanopore_only_genes), "\n")
cat("  Illumina-specific:", nrow(illumina_only_genes), "\n\n")

cat("Cross-Technology Correlation:\n")
cat("  Matched samples (median):", round(median(matched_correlations), 3), "\n")
cat("  All pairs (median):", round(median(cross_tech_corr), 3), "\n\n")

cat("Within-Technology Correlation:\n")
cat("  Nanopore (median):",
    round(median(nanopore_self_corr[lower.tri(nanopore_self_corr)]), 3), "\n")
cat("  Illumina (median):",
    round(median(illumina_self_corr[lower.tri(illumina_self_corr)]), 3), "\n\n")

# Save summary to file
summary_df <- data.frame(
  Metric = c(
    "Total_genes_nanopore",
    "Total_genes_illumina",
    "Overlapping_genes",
    "Nanopore_specific_genes",
    "Illumina_specific_genes",
    "Correlation_matched_median",
    "Correlation_all_pairs_median",
    "Correlation_nanopore_within_median",
    "Correlation_illumina_within_median",
    "Expr_pvalue_illumina",
    "Expr_pvalue_nanopore",
    "Length_pvalue_illumina",
    "Length_pvalue_nanopore"
  ),
  Value = c(
    nrow(nanopore_bed),
    nrow(illumina_bed),
    length(overlapping_genes),
    nrow(nanopore_only_genes),
    nrow(illumina_only_genes),
    median(matched_correlations),
    median(cross_tech_corr),
    median(nanopore_self_corr[lower.tri(nanopore_self_corr)]),
    median(illumina_self_corr[lower.tri(illumina_self_corr)]),
    pval_illumina,
    pval_nanopore,
    pval_length_illumina,
    pval_length_nanopore
  )
)

write.table(summary_df,
           "comparison_summary_statistics.txt",
           sep = "\t",
           quote = FALSE,
           row.names = FALSE)

cat("Analysis complete!\n")
cat("\nOutput files:\n")
cat("  - correlation_heatmap_nanopore_vs_illumina.pdf\n")
cat("  - matched_sample_correlations.pdf\n")
cat("  - all_pairs_correlations.pdf\n")
cat("  - nanopore_within_technology_correlation.pdf\n")
cat("  - illumina_within_technology_correlation.pdf\n")
cat("  - gene_expression_comparison_violin.pdf\n")
cat("  - gene_length_comparison_violin.pdf\n")
cat("  - comparison_summary_statistics.txt\n")
