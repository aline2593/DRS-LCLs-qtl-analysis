#!/usr/bin/env Rscript
################################################################################
# Complete FLAIR Downstream Analysis with SQANTI Reannotation
################################################################################
# Description: Comprehensive analysis from FLAIR quantification to final results:
#              1. Load and filter FLAIR isoforms
#              2. Run SQANTI3 for structural classification (via system call)
#              3. Reannotate using SQANTI results
#              4. Filter to protein-coding/lincRNA at 50% expression
#              5. Generate final publication plots (figure 1E and S2A-D, S3B)
#
# Input:  - FLAIR TPM quantification (per chromosome)
#         - Gencode annotation (GTF)
#         - Gene-level expression data
#
# Output: - Filtered BED file for tQTL analysis
#         - Reannotated isoform tables
#         - Publication-ready plots (figure 1E and S2A-D, S3B)
#
# Dependencies: data.table, parallel, dplyr, readr, stringr, ggplot2, tidyr
################################################################################

library(data.table)
library(parallel)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(gplots)

# Set working directory
workDir <- "/path/to/flair/output"
setwd(workDir)

cat("\n", rep("=", 70), "\n", sep = "")
cat("FLAIR COMPLETE DOWNSTREAM ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

################################################################################
# PART 1: LOAD AND FILTER FLAIR QUANTIFICATION
################################################################################

cat("PART 1: Loading and filtering FLAIR quantification...\n")
cat(rep("-", 70), "\n", sep = "")

################################################################################
# STEP 1: LOAD GTF ANNOTATION (ONCE, AT THE BEGINNING)
################################################################################

cat("Loading Gencode annotation...\n")
gtf_file <- "/path/to/annotation/gencode.v46lift37.annotation.gtf"
gtf_data <- fread(gtf_file, header = FALSE, sep = "\t")

# Extract gene-level information
genes_annotation <- gtf_data[V3 == "gene", .(
  chr = V1,
  start = V4,
  end = V5,
  gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
  gene_type = sub('.*gene_type "([^"]+)".*', '\\1', V9)
)]

cat("Total genes in annotation:", nrow(genes_annotation), "\n\n")

################################################################################
# STEP 2: LOAD AND MERGE CHROMOSOME-LEVEL TPM FILES
################################################################################

cat("Loading FLAIR TPM quantification files...\n")

# Define chromosomes to process
chromosomes <- c(1:22, "X", "Y", "M")

# Read TPM files for each chromosome
tpm_list <- lapply(chromosomes, function(chr) {
  file_path <- file.path(workDir, paste0("chr", chr, ".tpm.tsv"))
  cat("  Reading:", basename(file_path), "\n")
  read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
})

# Combine all chromosomes into single table
tpm_combined <- do.call(rbind, tpm_list)

# Clean column names (remove batch/condition labels if present)
colnames(tpm_combined) <- gsub("_condition.*|_batch.*", "", colnames(tpm_combined))

cat("\nTotal isoforms across all chromosomes:", nrow(tpm_combined), "\n")
cat("Number of samples:", ncol(tpm_combined) - 1, "\n")

# Check for duplicates
n_duplicates <- sum(duplicated(tpm_combined$ids))
if (n_duplicates > 0) {
  cat("WARNING: Found", n_duplicates, "duplicate isoforms - removing...\n")
  tpm_combined <- tpm_combined[!duplicated(tpm_combined), ]
}

# Save unfiltered data
write.table(tpm_combined,
            "01_all_isoforms_unfiltered.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Saved: 01_all_isoforms_unfiltered.txt\n\n")

################################################################################
# STEP 3: CLASSIFY ISOFORMS AS NOVEL VS ANNOTATED (INITIAL)
################################################################################

cat("Classifying isoforms as novel vs annotated...\n")

# Identify unknown/unannotated transcripts (FLAIR labels with "_chr")
unknown_idx <- grep("_chr", tpm_combined$ids)
cat("  Unknown/unannotated isoforms:", length(unknown_idx), "\n")

# Remove unknown isoforms
tpm_filtered <- tpm_combined[-unknown_idx, ]
cat("  Isoforms after removing unknowns:", nrow(tpm_filtered), "\n")

# Classify as annotated (ENST ID) or novel
tpm_filtered$classification <- ifelse(
  grepl("ENST", tpm_filtered$ids),
  "annotated",
  "novel"
)

cat("  Annotated isoforms:", sum(tpm_filtered$classification == "annotated"), "\n")
cat("  Novel isoforms:", sum(tpm_filtered$classification == "novel"), "\n\n")

# Save filtered data
write.table(tpm_filtered,
            "02_isoforms_no_unknowns.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# STEP 4: FILTER BY EXPRESSION THRESHOLD
################################################################################

cat("Filtering by expression threshold (TPM >= 5 in at least one sample)...\n")

# Count samples with TPM >= 5 for each isoform
n_samples_above_threshold <- mclapply(1:nrow(tpm_filtered), function(i) {
  sum(tpm_filtered[i, -c(1, ncol(tpm_filtered))] >= 5)
}, mc.cores = 6)

# Create summary table
isoform_summary <- data.frame(
  isoform_id = tpm_filtered$ids,
  n_samples_5tpm = unlist(n_samples_above_threshold),
  classification = tpm_filtered$classification,
  stringsAsFactors = FALSE
)

# Remove isoforms with TPM < 5 in all samples
isoform_summary_filtered <- isoform_summary[isoform_summary$n_samples_5tpm > 0, ]

cat("Isoforms with TPM >= 5 in at least one sample:", nrow(isoform_summary_filtered), "\n")
cat("  Annotated:", sum(isoform_summary_filtered$classification == "annotated"), "\n")
cat("  Novel:", sum(isoform_summary_filtered$classification == "novel"), "\n\n")

# Save summary
write.table(isoform_summary_filtered,
            "03_isoform_expression_summary.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# STEP 5: COUNT ISOFORMS PER GENE (BEFORE SQANTI)
################################################################################

cat("Counting isoforms per gene (before SQANTI reannotation)...\n")

# Standardize gene IDs (replace underscores with hyphens for matching)
genes_annotation$gene_id_std <- gsub("_", "-", genes_annotation$gene_id)

# Count annotated and novel isoforms per gene
cat("Processing genes (this may take a few minutes)...\n")

isoform_counts_per_gene <- mclapply(genes_annotation$gene_id_std, function(gene) {
  # Find isoforms matching this gene
  matching_isoforms <- grep(gene, isoform_summary_filtered$isoform_id)

  if (length(matching_isoforms) == 0) {
    return(list(n_annotated = 0, n_novel = 0))
  }

  # Count annotated (ENST) vs novel
  annotated_idx <- grep("ENST", isoform_summary_filtered$isoform_id[matching_isoforms])
  n_annotated <- length(annotated_idx)
  n_novel <- length(matching_isoforms) - n_annotated

  return(list(n_annotated = n_annotated, n_novel = n_novel))
}, mc.cores = 6)

# Compile results
gene_isoform_counts_before_sqanti <- data.frame(
  gene_id = genes_annotation$gene_id,
  n_annotated = sapply(isoform_counts_per_gene, function(x) x$n_annotated),
  n_novel = sapply(isoform_counts_per_gene, function(x) x$n_novel),
  stringsAsFactors = FALSE
)

# Filter to genes with at least one isoform
gene_isoform_counts_expressed <- gene_isoform_counts_before_sqanti[
  rowSums(gene_isoform_counts_before_sqanti[, c("n_annotated", "n_novel")]) > 0,
]

cat("Genes with at least one expressed isoform:", nrow(gene_isoform_counts_expressed), "\n\n")

# Save results
write.table(gene_isoform_counts_expressed,
            "04_isoforms_per_gene_before_sqanti.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# PART 2: PREPARE GTF FOR SQANTI AND RUN SQANTI3
################################################################################

cat("\nPART 2: Preparing files and running SQANTI3...\n")
cat(rep("-", 70), "\n", sep = "")

cat("Creating filtered GTF for SQANTI...\n")

# Extract transcript-level entries
transcripts_gtf <- gtf_data[V3 == "transcript"]

# Filter to protein-coding and lincRNA only
transcripts_filtered <- transcripts_gtf[
  grepl('gene_type "protein_coding"|gene_type "lincRNA"', V9)
]

# Extract transcript IDs that passed expression filter
filtered_transcript_ids <- gsub("_.*", "", isoform_summary_filtered$isoform_id)

# Keep only GTF entries for filtered transcripts
gtf_for_sqanti <- gtf_data[
  grepl(paste(filtered_transcript_ids, collapse = "|"), V9)
]

# Save filtered GTF
fwrite(gtf_for_sqanti,
       "Flair_filtered_for_sqanti.gtf",
       sep = "\t", col.names = FALSE, quote = FALSE)

cat("Filtered GTF created: Flair_filtered_for_sqanti.gtf\n")
cat("GTF entries:", nrow(gtf_for_sqanti), "\n\n")

# Run SQANTI3
cat("Running SQANTI3 classification...\n")
cat("This may take 1-3 hours...\n\n")

sqanti_cmd <- sprintf(
  "ml sqanti3 && python /path/to/SQANTI3/sqanti3_qc.py \\
    Flair_filtered_for_sqanti.gtf \\
    %s \\
    /path/to/reference/hg19.fa \\
    --force_id_ignore \\
    --polyA_motif_list /path/to/polyA_motifs/mouse_and_human.polyA_motif.txt \\
    -d sqanti_output \\
    --report both",
  gtf_file
)

# Execute SQANTI3
system(sqanti_cmd)

cat("SQANTI3 complete!\n\n")

# Check if SQANTI output exists
sqanti_classification_file <- "sqanti_output/Flair_filtered_for_sqanti_classification.txt"
if (!file.exists(sqanti_classification_file)) {
  stop("SQANTI classification file not found. Check SQANTI3 output.")
}

################################################################################
# PART 3: LOAD SQANTI RESULTS AND REANNOTATE
################################################################################

cat("PART 3: Reannotating isoforms using SQANTI classification...\n")
cat(rep("-", 70), "\n", sep = "")

# Load SQANTI classification
sqanti_class <- read_tsv(sqanti_classification_file,
                         comment = "#",
                         col_types = cols())

cat("SQANTI classifications loaded:", nrow(sqanti_class), "\n")

# Standardize isoform IDs
sqanti_class <- sqanti_class %>%
  mutate(
    isoform = if_else(
      str_starts(isoform, "ENST"),
      str_extract(isoform, "^[^_]+"),
      isoform
    )
  )

isoform_summary_filtered <- isoform_summary_filtered %>%
  mutate(isoform = str_extract(isoform_id, "^[^_]+"))

# Remove duplicates from SQANTI
sqanti_unique <- sqanti_class %>%
  distinct(isoform, .keep_all = TRUE)

# Merge FLAIR with SQANTI
merged_data <- isoform_summary_filtered %>%
  left_join(
    sqanti_unique %>%
      select(isoform, structural_category, associated_gene, associated_transcript),
    by = "isoform"
  )

cat("Merged FLAIR + SQANTI:", nrow(merged_data), "\n")

# Show initial classification
cat("\nBefore reannotation:\n")
print(table(merged_data$classification))

# Reannotate: novel isoforms with ENST associated_transcript → annotated
merged_data <- merged_data %>%
  mutate(
    classification_refined = case_when(
      classification == "novel" & str_detect(associated_transcript, "^ENST") ~ "annotated",
      TRUE ~ classification
    )
  )

cat("\nAfter reannotation:\n")
print(table(merged_data$classification_refined))

# Show reannotation statistics
reannotated_counts <- merged_data %>%
  filter(classification == "novel" & classification_refined == "annotated") %>%
  count(structural_category, name = "count") %>%
  arrange(desc(count))

cat("\nReannotated isoforms by structural category:\n")
print(reannotated_counts)

################################################################################
# PART 4: REMOVE DUPLICATE TRANSCRIPT IDs
################################################################################

cat("\n")
cat("PART 4: Removing duplicate transcript IDs...\n")
cat(rep("-", 70), "\n", sep = "")

# Find duplicate ENST IDs
annotated_enst <- merged_data %>%
  filter(classification_refined == "annotated",
         grepl("^ENST", associated_transcript)) %>%
  pull(associated_transcript) %>%
  unique()

reannotated_enst <- merged_data %>%
  filter(classification == "novel",
         classification_refined == "annotated",
         grepl("^ENST", associated_transcript)) %>%
  pull(associated_transcript) %>%
  unique()

duplicate_enst <- intersect(annotated_enst, reannotated_enst)

cat("Duplicate ENST IDs:", length(duplicate_enst), "\n")

# Remove duplicates (keep original annotated)
merged_data_clean <- merged_data %>%
  filter(!(associated_transcript %in% duplicate_enst &
           classification == "novel" &
           classification_refined == "annotated"))

cat("Isoforms after removing duplicates:", nrow(merged_data_clean), "\n")
cat("\nFinal classification:\n")
print(table(merged_data_clean$classification_refined))

# Save reannotated table
write.table(
  merged_data_clean %>% select(-classification),
  "05_isoforms_reannotated_sqanti.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

################################################################################
# PART 5: FILTER TO PROTEIN-CODING AND LINCRNA
################################################################################

cat("\n")
cat("PART 5: Filtering to protein-coding and lincRNA genes...\n")
cat(rep("-", 70), "\n", sep = "")

# Merge with gene types (using genes_annotation loaded at the beginning)
merged_with_type <- merged_data_clean %>%
  left_join(genes_annotation, by = c("associated_gene" = "gene_id"))

# Filter to protein-coding and lincRNA
protein_lincRNA <- merged_with_type %>%
  filter(gene_type %in% c("protein_coding", "lincRNA"))

cat("Protein-coding/lincRNA isoforms:", nrow(protein_lincRNA), "\n")

# Split by classification
novel_isoforms <- protein_lincRNA %>%
  filter(classification_refined == "novel")

annotated_isoforms <- protein_lincRNA %>%
  filter(classification_refined == "annotated")

cat("  Novel:", nrow(novel_isoforms), "\n")
cat("  Annotated:", nrow(annotated_isoforms), "\n\n")

# Save
write.table(novel_isoforms, "06_Novel_Proteincoding_lincRNA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(annotated_isoforms, "07_Annotated_Proteincoding_lincRNA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

################################################################################
# PART 6: PREPARE BED FILE FOR tQTL ANALYSIS (50% FILTER)
################################################################################

cat("PART 6: Preparing BED file for tQTL analysis...\n")
cat(rep("-", 70), "\n", sep = "")

# Extract transcript-level annotation
transcripts_annot <- gtf_data[V3 == "transcript", .(
  chr = V1,
  start = V4,
  end = V5,
  strand = V7,
  gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
  gene_type = sub('.*gene_type "([^"]+)".*', '\\1', V9),
  transcript_id = sub('.*transcript_id "([^"]+)".*', '\\1', V9)
)]

# Parse FLAIR IDs
isoforms_parsed <- tpm_filtered %>%
  separate(ids,
           into = c("transcript_id", "gene_id"),
           sep = "(?=ENSG)",
           extra = "merge",
           fill = "right")

isoforms_parsed$transcript_id <- sub("_([^_]*)$", "\\1", isoforms_parsed$transcript_id)
colnames(isoforms_parsed)[1:2] <- c("isoform_id", "gene_id")
isoforms_parsed$gene_id <- gsub("-", "_", isoforms_parsed$gene_id)

# Filter to protein-coding/lincRNA
transcripts_pc_lnc <- transcripts_annot[gene_type %in% c("protein_coding", "lincRNA")]

# Merge
transcripts_with_quant <- merge(
  transcripts_pc_lnc,
  isoforms_parsed,
  by = "gene_id",
  allow.cartesian = TRUE
)

# Remove duplicates
transcripts_unique <- transcripts_with_quant[!duplicated(transcripts_with_quant$isoform_id), ]

# Format as BED
bed_cols <- c("chr", "start", "end", "isoform_id", "gene_id", "strand")
sample_cols <- setdiff(colnames(transcripts_unique),
                       c(bed_cols, "gene_type", "transcript_id", "classification"))

transcripts_bed <- transcripts_unique[, c(bed_cols, sample_cols), with = FALSE]
colnames(transcripts_bed)[1:6] <- c("#chr", "start", "end", "id", "gid", "strd")

# Apply 50% expression filter
anno <- transcripts_bed[, 1:6]
quan <- transcripts_bed[, 7:ncol(transcripts_bed)]

zero_prop <- apply(quan, 1, function(x) sum(x == 0) / length(x))
keep <- zero_prop <= 0.5

transcripts_filtered <- cbind(anno[keep, ], quan[keep, ])

cat("Transcripts passing 50% filter:", nrow(transcripts_filtered), "\n")
cat("Unique genes:", length(unique(transcripts_filtered$gid)), "\n\n")

# Sort by position
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
transcripts_filtered$`#chr` <- factor(transcripts_filtered$`#chr`, levels = chr_order)

transcripts_sorted <- transcripts_filtered[
  order(transcripts_filtered$`#chr`,
        transcripts_filtered$start,
        transcripts_filtered$end),
]

# Save BED file for tQTL
write.table(transcripts_sorted,
            "08_transcripts_filtered50_for_tQTL.bed",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cat("BED file saved for tQTL analysis\n\n")

################################################################################
# PART 7: CALCULATE EXPRESSION STATISTICS
################################################################################

cat("PART 7: Calculating expression statistics...\n")
cat(rep("-", 70), "\n", sep = "")

# Calculate mean expression per transcript
tpm_with_class <- tpm_filtered %>%
  left_join(
    merged_data_clean %>% select(isoform_id, classification_refined, associated_gene),
    by = c("ids" = "isoform_id")
  )

# Calculate mean TPM
sample_columns <- setdiff(names(tpm_with_class),
                          c("ids", "classification", "classification_refined", "associated_gene"))

tpm_with_class$mean_tpm <- rowMeans(tpm_with_class[, sample_columns], na.rm = TRUE)
tpm_with_class$n_samples_5tpm <- apply(
  tpm_with_class[, sample_columns], 1,
  function(x) sum(x >= 5)
)

# Filter to protein-coding/lincRNA
novel_with_stats <- tpm_with_class %>%
  filter(classification_refined == "novel") %>%
  inner_join(genes_annotation, by = c("associated_gene" = "gene_id")) %>%
  filter(gene_type %in% c("protein_coding", "lincRNA"))

annotated_with_stats <- tpm_with_class %>%
  filter(classification_refined == "annotated") %>%
  inner_join(genes_annotation, by = c("associated_gene" = "gene_id")) %>%
  filter(gene_type %in% c("protein_coding", "lincRNA"))

cat("Expression statistics calculated\n\n")

################################################################################
# PART 8: GENERATE PUBLICATION PLOTS
################################################################################

cat("PART 8: Generating publication plots...\n")
cat(rep("-", 70), "\n", sep = "")

# Plot 1: Heatmap of annotated vs novel isoforms per gene
cat("Creating heatmap of isoforms per gene...\n")

# Count isoforms per gene (for protein-coding/lincRNA only, AFTER reannotation)
gene_isoform_counts_final <- protein_lincRNA %>%
  group_by(associated_gene, classification_refined) %>%
  summarize(n_isoforms = n(), .groups = "drop") %>%
  pivot_wider(names_from = classification_refined,
              values_from = n_isoforms,
              values_fill = 0) %>%
  rename(gene_id = associated_gene,
         n_annotated = annotated,
         n_novel = novel)

# Create contingency table
isoform_table <- table(
  gene_isoform_counts_final$n_annotated,
  gene_isoform_counts_final$n_novel
)

# Color palette and breaks
color_palette <- colorRampPalette(c("white", "orange", "red"))(n = 299)
color_breaks <- c(
  seq(0, 10, length = 100),
  seq(11, 100, length = 100),
  seq(101, 3000, length = 100)
)

# Create heatmap
pdf("09_heatmap_annotated_vs_novel_per_gene_1E.pdf", width = 10, height = 8)
heatmap.2(
  isoform_table,
  dendrogram = 'none',
  Rowv = FALSE,
  Colv = FALSE,
  trace = 'none',
  col = color_palette,
  breaks = color_breaks,
  key = FALSE,
  cellnote = ifelse(isoform_table == 0, NA, isoform_table),
  notecol = "black",
  xlab = "Number of Novel Isoforms",
  ylab = "Number of Annotated Isoforms",
  main = "Isoforms per Gene (Protein-coding & lincRNA, 50% filter)",
  cexRow = 1.2,
  cexCol = 1.2
)
dev.off()

# Plot 2: Transcript expression distribution
cat("Creating transcript expression plot...\n")

p_expr <- ggplot() +
  geom_histogram(aes(x = novel_with_stats$mean_tpm),
                 color = "grey28", fill = "cyan4", alpha = 1, bins = 50) +
  geom_histogram(aes(x = annotated_with_stats$mean_tpm),
                 color = "grey28", fill = "orange1", alpha = 0.7, bins = 50) +
  geom_vline(aes(xintercept = mean(novel_with_stats$mean_tpm, na.rm = TRUE)),
             color = "darkcyan", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean(annotated_with_stats$mean_tpm, na.rm = TRUE)),
             color = "orange3", linetype = "dashed", linewidth = 1) +
  scale_x_log10() +
  theme_light(base_size = 18) +
  labs(
    title = "Transcript Expression: Novel vs Annotated",
    x = "Mean TPM (log10)",
    y = "Frequency"
  ) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave("10_transcript_expression_distribution_S2A.pdf", plot = p_expr, height = 8, width = 10)

# Plot 3: Sample breadth
cat("Creating sample breadth plot...\n")

p_breadth <- ggplot() +
  geom_histogram(aes(x = novel_with_stats$n_samples_5tpm),
                 color = "grey28", fill = "cyan4", alpha = 1, bins = 30) +
  geom_histogram(aes(x = annotated_with_stats$n_samples_5tpm),
                 color = "grey28", fill = "orange1", alpha = 0.7, bins = 30) +
  geom_vline(aes(xintercept = mean(novel_with_stats$n_samples_5tpm, na.rm = TRUE)),
             color = "darkcyan", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean(annotated_with_stats$n_samples_5tpm, na.rm = TRUE)),
             color = "orange3", linetype = "dashed", linewidth = 1) +
  theme_light(base_size = 18) +
  labs(
    title = "Sample Breadth: Novel vs Annotated Transcripts",
    x = "Number of Samples with TPM >= 5",
    y = "Frequency"
  ) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave("11_transcript_sample_breadth_S3B.pdf", plot = p_breadth, height = 8, width = 10)

# Plot 4-6: Gene expression vs isoform number
cat("Creating gene expression vs isoform plots...\n")

# Load gene expression
gene_expr <- fread("/path/to/gene/expression/rpkm_edgeR.txt")
gene_expr_mean <- data.frame(
  gene_id = gene_expr[[1]],
  mean_gene_expr = rowMeans(gene_expr[, -1])
)

# Merge with expression
gene_analysis <- gene_isoform_counts_final %>%
  left_join(gene_expr_mean, by = "gene_id") %>%
  filter(!is.na(mean_gene_expr)) %>%
  mutate(total_isoforms = n_annotated + n_novel)

cat("Genes with expression and isoform data:", nrow(gene_analysis), "\n")

# Plot: Annotated isoforms
p_anno <- ggplot(gene_analysis, aes(x = as.factor(n_annotated), y = mean_gene_expr)) +
  geom_violin(fill = "orange1", trim = FALSE) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.3) +
  scale_y_log10() +
  theme_light(base_size = 14) +
  labs(
    x = "Number of Annotated Isoforms",
    y = "Mean Gene Expression (log10 TPM)"
  ) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave("12_gene_expr_vs_annotated_isoforms_S2B.pdf", plot = p_anno, height = 10, width = 20)

# Plot: Novel isoforms
p_nov <- ggplot(gene_analysis, aes(x = as.factor(n_novel), y = mean_gene_expr)) +
  geom_violin(fill = "cyan4", trim = FALSE) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.3) +
  scale_y_log10() +
  theme_light(base_size = 14) +
  labs(
    x = "Number of Novel Isoforms",
    y = "Mean Gene Expression (log10 TPM)"
  ) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave("13_gene_expr_vs_novel_isoforms_S2C.pdf", plot = p_nov, height = 10, width = 20)

# Plot: Total isoforms
p_tot <- ggplot(gene_analysis, aes(x = as.factor(total_isoforms), y = mean_gene_expr)) +
  geom_violin(fill = "tomato3", trim = FALSE) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.3) +
  scale_y_log10() +
  theme_light(base_size = 14) +
  labs(
    x = "Total Number of Isoforms",
    y = "Mean Gene Expression (log10 TPM)"
  ) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave("14_gene_expr_vs_total_isoforms_S2D.pdf", plot = p_tot, height = 10, width = 25)

cat("All plots generated!\n\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat(rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Summary Statistics:\n")
cat("  Initial isoforms (all chromosomes):", nrow(tpm_combined), "\n")
cat("  After removing unknowns:", nrow(tpm_filtered), "\n")
cat("  After TPM >= 5 filter:", nrow(isoform_summary_filtered), "\n")
cat("  After SQANTI reannotation:", nrow(merged_data_clean), "\n")
cat("  Duplicates removed:", length(duplicate_enst), "\n")
cat("  Final protein-coding/lincRNA:\n")
cat("    Novel:", nrow(novel_isoforms), "\n")
cat("    Annotated:", nrow(annotated_isoforms), "\n")
cat("  For tQTL (50% filter):", nrow(transcripts_sorted), "\n")
cat("  Final genes:", length(unique(transcripts_sorted$gid)), "\n\n")

cat("Expression Statistics:\n")
cat("  Novel - Mean TPM:", round(mean(novel_with_stats$mean_tpm, na.rm = TRUE), 2), "\n")
cat("  Annotated - Mean TPM:", round(mean(annotated_with_stats$mean_tpm, na.rm = TRUE), 2), "\n")
cat("  Novel - Mean samples >= 5 TPM:",
    round(mean(novel_with_stats$n_samples_5tpm, na.rm = TRUE), 2), "\n")
cat("  Annotated - Mean samples >= 5 TPM:",
    round(mean(annotated_with_stats$n_samples_5tpm, na.rm = TRUE), 2), "\n\n")

cat("Output Files:\n")
cat("  01_all_isoforms_unfiltered.txt\n")
cat("  02_isoforms_no_unknowns.txt\n")
cat("  03_isoform_expression_summary.txt\n")
cat("  04_isoforms_per_gene_before_sqanti.txt\n")
cat("  05_isoforms_reannotated_sqanti.txt\n")
cat("  06_Novel_Proteincoding_lincRNA.txt\n")
cat("  07_Annotated_Proteincoding_lincRNA.txt\n")
cat("  08_transcripts_filtered50_for_tQTL.bed  => USE FOR tQTL MAPPING\n")
cat("  09_heatmap_annotated_vs_novel_per_gene_1E.pdf\n")
cat("  10_transcript_expression_distribution_S2A.pdf\n")
cat("  11_transcript_sample_breadth_S3B.pdf\n")
cat("  12_gene_expr_vs_annotated_isoforms_S2B.pdf\n")
cat("  13_gene_expr_vs_novel_isoforms_S2C.pdf\n")
cat("  14_gene_expr_vs_total_isoforms_S2D.pdf\n")
cat("  sqanti_output/ (SQANTI3 results)\n\n")

cat("Next Step:\n")
cat("  bgzip 08_transcripts_filtered50_for_tQTL.bed\n")
cat("  tabix -p bed 08_transcripts_filtered50_for_tQTL.bed.gz\n")
cat(rep("=", 70), "\n", sep = "")
