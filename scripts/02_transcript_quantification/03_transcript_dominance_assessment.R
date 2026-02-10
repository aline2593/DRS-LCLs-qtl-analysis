#!/usr/bin/env Rscript
################################################################################
# Transcript Dominance Assessment
################################################################################
# Description: Determine if genes have a dominant transcript (one isoform
#              accounting for >90% of the gene's total expression) or exhibit
#              multi-isoform expression patterns.
#
# Methodology:
#   For each gene:
#   1. Calculate each transcript's percentage of total gene expression
#   2. Classify as "Dominant" if max transcript >90% of total
#   3. Classify as "No_dominant" if max transcript ≤90%
#
# Input:  - Protein-coding/lincRNA isoforms with 5TPM in at least one sample (no 50% filters)
#         - Gene-level expression counts
#
# Output: - Dominance classification per gene
#         - Comparative plots (dominant vs non-dominant genes)
#         - Statistics on isoform complexity
#
# Dependencies: data.table, parallel, dplyr, ggplot2, tibble
################################################################################

library(data.table)
library(parallel)
library(dplyr)
library(ggplot2)
library(tibble)

# Set working directory
workDir <- "/path/to/flair/output"
setwd(workDir)

cat("\n", rep("=", 70), "\n", sep = "")
cat("TRANSCRIPT DOMINANCE ASSESSMENT\n")
cat(rep("=", 70), "\n\n", sep = "")

################################################################################
# PART 1: LOAD DATA AND CALCULATE AVERAGE EXPRESSION
################################################################################

cat("PART 1: Loading isoform and gene expression data...\n")
cat(rep("-", 70), "\n", sep = "")

# Load protein-coding/lincRNA isoforms (from previous analysis)
# This should be the full set before 50% filter, with all samples
isoforms <- fread("Isoforms_protein_lincRNA_All.txt", data.table = FALSE)

cat("Total isoforms loaded:", nrow(isoforms), "\n")
cat("Total genes:", length(unique(isoforms$gid)), "\n\n")

# Extract TPM values (excluding coordinate columns)
# Columns: #chr, start, end, id (transcript), gid (gene), strd, sample1, sample2, ...
isoform_tpm <- isoforms[, 7:ncol(isoforms)]
rownames(isoform_tpm) <- isoforms$id

# Calculate mean TPM per transcript
mean_tpm_per_transcript <- data.frame(
  Transcript = isoforms$id,
  Gene = isoforms$gid,
  Average_transcripts_tpm = rowMeans(isoform_tpm)
)

cat("Calculated mean TPM per transcript\n")

################################################################################
# PART 2: LOAD TRANSCRIPT-LEVEL COUNTS
################################################################################

cat("\nPART 2: Loading transcript-level raw counts...\n")
cat(rep("-", 70), "\n", sep = "")

# Load FLAIR raw counts (not TPM)
count_transcripts <- fread("/path/to/flair/counts/all_chr.counts.tsv",
                           data.table = FALSE)

cat("Total transcript count entries:", nrow(count_transcripts), "\n")

# Separate annotated and novel transcripts

# Annotated (contain ENST)
annotated_idx <- grep("ENST", count_transcripts$ids)
annotated_counts <- count_transcripts[annotated_idx, ]

cat("Annotated transcripts:", nrow(annotated_counts), "\n")

# Calculate mean counts for annotated
annotated_counts_only <- annotated_counts[, -1]
rownames(annotated_counts_only) <- annotated_counts$ids

annotated_mean <- data.frame(
  ids = annotated_counts$ids,
  Average_transcripts_count = rowMeans(annotated_counts_only)
)

# Parse transcript and gene IDs from annotated (format: ENST_ENSG)
annotated_parsed <- annotated_mean %>%
  separate(ids, into = c("Transcript", "Gene"), sep = "_ENSG", remove = FALSE) %>%
  mutate(Gene = paste0("ENSG", Gene)) %>%
  select(Transcript, Gene, Average_transcripts_count)

# Novel (do not contain ENST)
novel_idx <- grep("ENST", count_transcripts$ids, invert = TRUE)
novel_counts <- count_transcripts[novel_idx, ]
novel_counts <- novel_counts[!duplicated(novel_counts$ids), ]

cat("Novel transcripts:", nrow(novel_counts), "\n")

# Calculate mean counts for novel
novel_counts_only <- novel_counts[, -1]
rownames(novel_counts_only) <- novel_counts$ids

novel_mean <- data.frame(
  ids = novel_counts$ids,
  Average_transcripts_count = rowMeans(novel_counts_only)
)

# Parse transcript and gene IDs from novel (format: transcript_ENSG)
# Keep only those with ENSG (filter out incomplete IDs)
novel_with_gene <- novel_mean[grep("ENSG", novel_mean$ids), ]

novel_parsed <- novel_with_gene %>%
  separate(ids, into = c("Transcript", "rest"), sep = "_", extra = "merge", remove = FALSE) %>%
  separate(rest, into = c("Gene", "extra"), sep = "_", extra = "merge", remove = TRUE) %>%
  select(Transcript, Gene, Average_transcripts_count)

cat("Novel transcripts with gene assignment:", nrow(novel_parsed), "\n\n")

################################################################################
# PART 3: LOAD GENE-LEVEL COUNTS
################################################################################

cat("PART 3: Loading gene-level counts...\n")
cat(rep("-", 70), "\n", sep = "")

# Load gene-level counts (from gene expression analysis)
gene_counts <- fread("/path/to/gene/expression/Count_all.txt.gz",
                    data.table = FALSE)

cat("Total genes with counts:", nrow(gene_counts), "\n")

# Calculate mean gene counts
gene_counts_only <- gene_counts[, -ncol(gene_counts)]  # Remove length column if present
gene_id_col <- gene_counts_only[, 1]
gene_counts_matrix <- gene_counts_only[, -1]

mean_gene_counts <- data.frame(
  Gene = gene_id_col,
  Average_gene_count = rowMeans(gene_counts_matrix)
)

cat("Calculated mean counts per gene\n\n")

################################################################################
# PART 4: MERGE TRANSCRIPT TPM AND COUNTS
################################################################################

cat("PART 4: Merging transcript data...\n")
cat(rep("-", 70), "\n", sep = "")

# Merge annotated transcripts
annotated_complete <- merge(
  mean_tpm_per_transcript,
  annotated_parsed,
  by = "Transcript"
) %>%
  select(Transcript, Gene = Gene.x, Average_transcripts_tpm, Average_transcripts_count)

cat("Annotated transcripts merged:", nrow(annotated_complete), "\n")

# Merge novel transcripts
novel_complete <- merge(
  mean_tpm_per_transcript,
  novel_parsed,
  by = "Transcript"
) %>%
  select(Transcript, Gene = Gene.x, Average_transcripts_tpm, Average_transcripts_count)

cat("Novel transcripts merged:", nrow(novel_complete), "\n")

# Add gene-level counts
annotated_with_gene <- merge(annotated_complete, mean_gene_counts, by = "Gene")
annotated_with_gene$group <- "Annotated"

novel_with_gene <- merge(novel_complete, mean_gene_counts, by = "Gene")
novel_with_gene$group <- "Novel"

# Combine
all_transcripts <- rbind(annotated_with_gene, novel_with_gene)

cat("Total transcripts with complete data:", nrow(all_transcripts), "\n")
cat("Total genes:", length(unique(all_transcripts$Gene)), "\n\n")

write.table(all_transcripts,
            "01_transcripts_with_counts_and_tpm.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# PART 5: CALCULATE DOMINANCE
################################################################################

cat("PART 5: Calculating transcript dominance...\n")
cat(rep("-", 70), "\n", sep = "")

unique_genes <- unique(all_transcripts$Gene)
cat("Processing", length(unique_genes), "genes...\n\n")

# For each gene, calculate percentage of total counts per transcript
transcript_percentages <- mclapply(unique_genes, function(gene) {
  # Get all transcripts for this gene
  gene_transcripts <- all_transcripts[all_transcripts$Gene == gene, ]

  # Calculate percentage each transcript contributes
  total_counts <- sum(gene_transcripts$Average_transcripts_count)

  if (total_counts == 0) {
    percentages <- rep(0, nrow(gene_transcripts))
  } else {
    percentages <- (gene_transcripts$Average_transcripts_count / total_counts) * 100
  }

  return(percentages)
}, mc.cores = 6)

# Find maximum percentage per gene
max_percentages <- sapply(transcript_percentages, max)

cat("Genes with single transcript (100%):", sum(max_percentages == 100), "\n")
cat("Genes with dominant transcript (>90% but <100%):",
    sum(max_percentages > 90 & max_percentages < 100), "\n")
cat("Genes with no dominant transcript (≤90%):", sum(max_percentages <= 90), "\n\n")

# Count transcripts per gene
n_transcripts_per_gene <- sapply(transcript_percentages, length)

# Find position of max transcript
max_transcript_idx <- sapply(transcript_percentages, which.max)

# Get transcript ID of dominant/highest abundance transcript
dominant_transcript_ids <- sapply(1:length(unique_genes), function(i) {
  gene_transcripts <- all_transcripts[all_transcripts$Gene == unique_genes[i], ]
  gene_transcripts$Transcript[max_transcript_idx[i]]
})

# Compile results
dominance_results <- data.frame(
  gene = unique_genes,
  maxTranscript = dominant_transcript_ids,
  maxTranscriptPercentage = max_percentages,
  numberOfTranscripts = n_transcripts_per_gene,
  stringsAsFactors = FALSE
)

cat("Dominance calculated for all genes\n\n")

################################################################################
# PART 6: CLASSIFY GENES AS DOMINANT OR NOT
################################################################################

cat("PART 6: Classifying genes...\n")
cat(rep("-", 70), "\n", sep = "")

# Classify: >90% = Dominant
dominant_genes <- dominance_results %>%
  filter(maxTranscriptPercentage > 90)

no_dominant_genes <- dominance_results %>%
  filter(maxTranscriptPercentage <= 90)

cat("Dominant genes (>90%):", nrow(dominant_genes), "\n")
cat("No dominant genes (≤90%):", nrow(no_dominant_genes), "\n\n")

# Save classification
write.table(dominant_genes,
            "02_genes_with_dominant_transcript.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(no_dominant_genes,
            "03_genes_without_dominant_transcript.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# PART 7: CONVERT GENE COUNTS TO TPM FOR COMPARISON
################################################################################

cat("PART 7: Converting gene-level counts to TPM...\n")
cat(rep("-", 70), "\n", sep = "")

# Load gene RPKM
gene_rpkm <- fread("/path/to/gene/expression/rpkm_edgeR.txt", data.table = FALSE)

# Convert RPKM to TPM (normalize each sample to sum to 1 million)
gene_ids <- gene_rpkm[[1]]
rpkm_matrix <- as.matrix(gene_rpkm[, -1])

tpm_matrix <- apply(rpkm_matrix, 2, function(x) {
  (x / sum(x, na.rm = TRUE)) * 1e6
})

gene_tpm <- data.table(Geneid = gene_ids, tpm_matrix)

# Calculate mean TPM per gene (exclude zero-only genes)
gene_tpm_matrix <- as.matrix(gene_tpm[, -1])
rownames(gene_tpm_matrix) <- gene_tpm$Geneid

# Remove genes with zero expression in all samples
keep_genes <- rowSums(gene_tpm_matrix > 0) > 0
gene_tpm_nonzero <- gene_tpm_matrix[keep_genes, ]

cat("Genes with expression:", nrow(gene_tpm_nonzero), "\n")

# Calculate mean
mean_gene_tpm <- data.frame(
  gene = rownames(gene_tpm_nonzero),
  Average_gene_tpm = rowMeans(gene_tpm_nonzero)
)

# Merge with dominance classification
dominant_with_expr <- merge(mean_gene_tpm, dominant_genes, by = "gene")
dominant_with_expr$status <- "Dominant"

no_dominant_with_expr <- merge(mean_gene_tpm, no_dominant_genes, by = "gene")
no_dominant_with_expr$status <- "No_dominant"

cat("Dominant genes with expression:", nrow(dominant_with_expr), "\n")
cat("No dominant genes with expression:", nrow(no_dominant_with_expr), "\n\n")

################################################################################
# PART 8: GENERATE PLOTS
################################################################################

cat("PART 8: Generating plots...\n")
cat(rep("-", 70), "\n", sep = "")

# Plot 1: Gene expression comparison (dominant vs no dominant)
cat("Creating gene expression comparison plot...\n")

p_expr <- ggplot() +
  geom_histogram(aes(x = dominant_with_expr$Average_gene_tpm),
                 color = "grey28", fill = "darkmagenta", alpha = 1, bins = 50) +
  geom_histogram(aes(x = no_dominant_with_expr$Average_gene_tpm),
                 color = "grey28", fill = "darkolivegreen", alpha = 0.7, bins = 50) +
  geom_vline(aes(xintercept = mean(dominant_with_expr$Average_gene_tpm)),
             color = "darkred", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean(no_dominant_with_expr$Average_gene_tpm)),
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  scale_x_log10() +
  theme_light(base_size = 18) +
  labs(
    title = "Gene Expression: Dominant vs No Dominant Transcript",
    x = "Mean Gene Expression (log10 TPM)",
    y = "Frequency"
  ) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  )

ggsave("04_gene_expression_dominant_vs_nodominant_S3C.pdf",
       plot = p_expr, height = 8, width = 10)

# Statistical test
wilcox_test <- wilcox.test(dominant_with_expr$Average_gene_tpm,
                           no_dominant_with_expr$Average_gene_tpm)

cat("Gene expression comparison:\n")
cat("  Dominant mean TPM:", round(mean(dominant_with_expr$Average_gene_tpm), 2), "\n")
cat("  No dominant mean TPM:", round(mean(no_dominant_with_expr$Average_gene_tpm), 2), "\n")
cat("  Wilcoxon p-value:", format(wilcox_test$p.value, scientific = TRUE), "\n\n")

################################################################################
# PART 9: CREATE FINAL COMBINED TABLE
################################################################################

cat("PART 9: Creating final combined table...\n")
cat(rep("-", 70), "\n", sep = "")

# Combine dominant and no dominant with all info
combined_table <- rbind(
  dominant_with_expr[, c("gene", "Average_gene_tpm", "maxTranscript",
                         "maxTranscriptPercentage", "numberOfTranscripts", "status")],
  no_dominant_with_expr[, c("gene", "Average_gene_tpm", "maxTranscript",
                           "maxTranscriptPercentage", "numberOfTranscripts", "status")]
)

cat("Combined table:", nrow(combined_table), "genes\n")

write.table(combined_table,
            "07_gene_dominance_classification.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("DOMINANCE ANALYSIS COMPLETE!\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Summary Statistics:\n")
cat("  Total genes analyzed:", length(unique_genes), "\n")
cat("  Genes with single transcript:", sum(max_percentages == 100), "\n")
cat("  Genes with dominant transcript (>90%):", nrow(dominant_genes), "\n")
cat("  Genes without dominant transcript (≤90%):", nrow(no_dominant_genes), "\n\n")

cat("Expression Characteristics:\n")
cat("  Dominant genes - Mean TPM:",
    round(mean(dominant_with_expr$Average_gene_tpm), 2), "\n")
cat("  No dominant genes - Mean TPM:",
    round(mean(no_dominant_with_expr$Average_gene_tpm), 2), "\n")
cat("  Wilcoxon test p-value:",
    format(wilcox_test$p.value, scientific = TRUE), "\n\n")

cat("Transcript Complexity:\n")
cat("  Dominant genes - Median transcripts:",
    median(dominant_genes$numberOfTranscripts), "\n")
cat("  No dominant genes - Median transcripts:",
    median(no_dominant_genes$numberOfTranscripts), "\n\n")

cat("Output Files:\n")
cat("  01_transcripts_with_counts_and_tpm.txt\n")
cat("  02_genes_with_dominant_transcript.txt\n")
cat("  03_genes_without_dominant_transcript.txt\n")
cat("  04_gene_expression_dominant_vs_nodominant.pdf\n")
cat("  05_dominant_transcript_by_complexity.pdf\n")
cat("  06_nodominant_transcript_by_complexity.pdf\n")
cat("  07_gene_dominance_classification.txt  =>  FINAL TABLE\n\n")

cat("Next Step:\n")
cat("  Use 07_gene_dominance_classification.txt for:\n")
cat("  - QTL stratification analysis\n")
cat("  - Functional enrichment\n")
cat("  - Correlation with gene features\n")
cat(rep("=", 70), "\n", sep = "")
