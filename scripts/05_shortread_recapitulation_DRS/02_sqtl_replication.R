#!/usr/bin/env Rscript
################################################################################
# sQTL Replication in Long-Read RNA-seq (Nanopore DRS)
################################################################################
# Description: Test replication of GTEx short-read sQTL discoveries in
#              long-read nanopore trQTL data from LCL samples
#
# Approach:
#   1. Lift over GTEx sQTL coordinates from hg38 to hg37 to match nanopore built
#   2. Match SNP positions to rsIDs using GTEx lookup table
#   3. Filter nanopore trQTL results to significant GTEx sQTL gene-SNP pairs
#   4. Calculate FDR for replication
#   5. Compare structural categories (SQANTI3) between overlapping/non-overlapping
#
# Input:  - GTEx v8 sQTL results (significant pairs + sgenes)
#         - GTEx variant lookup table (position → rsID mapping)
#         - Nanopore trQTL nominal results (all gene-transcript-SNP associations)
#         - SQANTI3 transcript classifications
#
# Output: - Replicated sQTLs in nanopore data at FDR < 5%
#         - Structural category comparison plots
#         - Summary statistics
#
# Dependencies: R (data.table, dplyr, ggplot2, qvalue)
################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(qvalue)

# Configuration
FDR_THRESHOLD <- 0.05

# Create output directories
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("======================================================================\n")
cat("sQTL Replication Analysis: GTEx (short-read) → Nanopore (long-read)\n")
cat("======================================================================\n\n")

################################################################################
# PART 1: COORDINATE LIFTOVER (hg38 → hg37)
################################################################################

cat("PART 1: Processing GTEx sQTL coordinates and rsID mapping\n")
cat("----------------------------------------------------------------------\n")

# Load lifted-over coordinates (hg38 → hg37)
# Generated using UCSC liftOver tool
sqtl_coords_lifted <- fread("hglft_v37_sQTL_SNP.txt")
cat("  ✓ Loaded lifted coordinates:", nrow(sqtl_coords_lifted), "variants\n")

# Load GTEx variant lookup table (position → rsID)
gtex_variant_lookup <- fread("GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt")
cat("  ✓ Loaded GTEx variant lookup:", nrow(gtex_variant_lookup), "variants\n")

# Clean chromosome names
sqtl_coords_lifted$V1 <- gsub("chr", "", sqtl_coords_lifted$V1)

# Merge to get rsIDs for lifted coordinates
sqtl_with_rsids <- merge(
  sqtl_coords_lifted,
  gtex_variant_lookup,
  by.x = c("V3", "V1"),  # v38 position, chromosome
  by.y = c("variant_pos", "chr")
)

cat("  ✓ Matched rsIDs:", nrow(sqtl_with_rsids), "variants\n")
cat("  ⚠ Lost", nrow(sqtl_coords_lifted) - nrow(sqtl_with_rsids),
    "variants in lookup\n")

# Save intermediate file
fwrite(sqtl_with_rsids,
       "results/01_sQTL_v37_v38_coordinates_with_rsIDs.txt",
       sep = "\t")

cat("\n")

################################################################################
# PART 2: LOAD GTEx sQTL RESULTS
################################################################################

cat("PART 2: Loading GTEx sQTL discovery results\n")
cat("----------------------------------------------------------------------\n")

# Load significant sQTL pairs from GTEx v8
gtex_sqtl_pairs <- fread("Cells_EBV-transformed_lymphocytes.v8.sqtl_signifpairs.txt")
cat("  ✓ Loaded significant sQTL pairs:", nrow(gtex_sqtl_pairs), "\n")
cat("  ✓ Unique phenotypes (introns):", uniqueN(gtex_sqtl_pairs$phenotype_id), "\n")

# Merge coordinates with sQTL results
gtex_sqtl_coords <- merge(
  gtex_sqtl_pairs,
  sqtl_with_rsids,
  by = c("chr", "position_v38")
)

# Remove duplicate positions (keep one per position)
gtex_sqtl_coords_unique <- gtex_sqtl_coords[!duplicated(position_v38)]

cat("  ✓ sQTL pairs with coordinates:", nrow(gtex_sqtl_coords_unique), "\n")

# Load sGene results (gene-level summary)
gtex_sgenes <- fread("Cells_EBV-transformed_lymphocytes.v8.sgenes_mod.txt")
cat("  ✓ Loaded sGenes:", nrow(gtex_sgenes), "\n")

# Filter to significant sGenes (FDR < 5%)
gtex_sgenes_sig <- gtex_sgenes %>%
  filter(pval_perm <= FDR_THRESHOLD)

cat("  ✓ Significant sGenes (permutation p ≤ 0.05):",
    nrow(gtex_sgenes_sig), "\n")

# Load phenotype groups (intron → gene mapping)
intron_to_gene <- fread(
  "Cells_EBV-transformed_lymphocytes.leafcutter.phenotype_groups.txt",
  header = FALSE
)
setnames(intron_to_gene, c("intron_id", "gene_id"))

cat("  ✓ Unique genes with introns:", uniqueN(intron_to_gene$gene_id), "\n")

# Merge everything together
gtex_sqtl_complete <- gtex_sqtl_coords_unique %>%
  left_join(intron_to_gene, by = c("phenotype_id" = "intron_id")) %>%
  left_join(
    gtex_sgenes %>% select(phenotype_id, gene_id, pval_perm, pval_nominal),
    by = c("phenotype_id", "gene_id")
  )

cat("  ✓ Complete sQTL dataset:", nrow(gtex_sqtl_complete), "associations\n")
cat("  ✓ Unique genes:", uniqueN(gtex_sqtl_complete$gene_id), "\n")

# Save complete sQTL set
fwrite(gtex_sqtl_complete,
       "results/02_GTEx_sQTL_complete_with_coordinates.txt",
       sep = "\t")

cat("\n")

################################################################################
# PART 3: PREPARE NANOPORE trQTL DATA
################################################################################

cat("PART 3: Loading nanopore trQTL results\n")
cat("----------------------------------------------------------------------\n")

# Load nanopore transcript quantifications
# Expected location: transcripts_annotated_novel_quant.filtered50.bed.gz
nanopore_transcripts <- fread("transcripts_annotated_novel_quant.filtered50.bed.gz")

# Clean gene IDs (remove version numbers)
nanopore_transcripts$gene_id <- sub("\\..*", "", nanopore_transcripts$gid)

cat("  ✓ Nanopore transcripts tested:", nrow(nanopore_transcripts), "\n")
cat("  ✓ Unique genes:", uniqueN(nanopore_transcripts$gene_id), "\n")

# Check overlap with GTEx sQTL genes
genes_gtex_sqtl <- unique(gtex_sqtl_complete$gene_id)
genes_nanopore <- unique(nanopore_transcripts$gene_id)
genes_common <- intersect(genes_gtex_sqtl, genes_nanopore)

cat("  ✓ sQTL genes tested in nanopore:", length(genes_common), "\n")
cat("  ⚠ sQTL genes NOT in nanopore:",
    length(genes_gtex_sqtl) - length(genes_common),
    "(not expressed at 50% threshold)\n\n")

# Filter GTEx sQTLs to genes tested in nanopore
gtex_sqtl_testable <- gtex_sqtl_complete %>%
  filter(gene_id %in% genes_common)

fwrite(gtex_sqtl_testable,
       "results/03_GTEx_sQTL_testable_in_nanopore.txt",
       sep = "\t")

cat("\n")

################################################################################
# PART 4: EXTRACT MATCHING trQTL RESULTS
################################################################################

cat("PART 4: Filtering nanopore trQTLs to GTEx sQTL gene-SNP pairs\n")
cat("----------------------------------------------------------------------\n")

# Create gene-rsID pairs from GTEx sQTLs
gtex_gene_snp_pairs <- gtex_sqtl_testable %>%
  select(gene_id, rs_id_dbSNP151_GRCh38p7) %>%
  distinct()

fwrite(gtex_gene_snp_pairs,
       "results/04_GTEx_gene_SNP_pairs_for_filtering.txt",
       sep = "\t",
       col.names = FALSE)

cat("  ✓ Created", nrow(gtex_gene_snp_pairs), "gene-SNP pairs for filtering\n")

# Note: The actual filtering should be done with bash/awk for memory efficiency
# Run the companion bash script: bash filter_nanopore_trqtls.sh
# Input: nominals_norm_all_combined.noheader.txt.gz (nanopore nominal results)
# Output: trQTL_nominals_matching_sQTL.txt.gz

if (!file.exists("trQTL_nominals_matching_sQTL.txt.gz")) {
  cat("\n  ⚠ ERROR: trQTL_nominals_matching_sQTL.txt.gz not found!\n")
  cat("  ⚠ Please run: bash filter_nanopore_trqtls.sh\n")
  cat("  ⚠ This script filters the large nominal results file\n\n")
  stop("Missing filtered trQTL file. Run bash filtering script first.")
}

cat("  → Loading pre-filtered nanopore results...\n")

# Load the filtered results
nanopore_trqtl_filtered <- fread("trQTL_nominals_matching_sQTL.txt.gz",
                                  header = FALSE)

# Set column names (QTLtools nominal output format)
setnames(nanopore_trqtl_filtered,
         c("grp_id", "phe_chr", "phe_from", "phe_to", "phe_strand",
           "phe_id", "n_phe_in_grp", "n_var_in_cis", "dist_phe_var",
           "var_id", "var_chr", "var_from", "var_to",
           "nom_pval", "r_squared", "slope", "best_hit")[1:ncol(nanopore_trqtl_filtered)])

cat("  ✓ Loaded filtered nanopore trQTLs:", nrow(nanopore_trqtl_filtered), "\n")
cat("  ✓ Unique genes:", uniqueN(nanopore_trqtl_filtered$grp_id), "\n")
cat("  ✓ Unique SNPs:", uniqueN(nanopore_trqtl_filtered$var_id), "\n")

# Clean gene IDs
nanopore_trqtl_filtered$gene_id <- sub("\\..*", "", nanopore_trqtl_filtered$grp_id)

cat("\n")

################################################################################
# PART 5: MERGE AND CALCULATE FDR
################################################################################

cat("PART 5: sQTL replication analysis in nanopore trQTLs\n")
cat("----------------------------------------------------------------------\n")

# Merge nanopore trQTLs with GTEx sQTLs
sqtl_trqtl_matched <- merge(
  nanopore_trqtl_filtered,
  gtex_sqtl_testable,
  by.x = c("gene_id", "var_id"),
  by.y = c("gene_id", "rs_id_dbSNP151_GRCh38p7")
)

# Remove duplicate phenotypes
sqtl_trqtl_unique <- sqtl_trqtl_matched[!duplicated(phenotype_id)]

cat("  ✓ Matched associations:", nrow(sqtl_trqtl_unique), "\n")
cat("  ✓ Unique genes:", uniqueN(sqtl_trqtl_unique$gene_id), "\n")
cat("  ✓ Unique GTEx introns:", uniqueN(sqtl_trqtl_unique$phenotype_id), "\n")
cat("  ✓ Unique nanopore transcripts:", uniqueN(sqtl_trqtl_unique$phe_id), "\n")
cat("  ✓ Unique SNPs:", uniqueN(sqtl_trqtl_unique$var_id), "\n\n")

# Calculate q-values
cat("  Calculating FDR (q-values)...\n")

qobj <- qvalue(p = sqtl_trqtl_unique$nom_pval, fdr.level = FDR_THRESHOLD)

cat("\n  === sQTL Replication Summary ===\n")
cat("  π₀ (proportion of true nulls):", round(qobj$pi0, 3), "\n")
cat("  π₁ (proportion of true signals):", round(1 - qobj$pi0, 3), "\n\n")

print(summary(qobj))

# Extract significant replications
sqtl_replicated <- sqtl_trqtl_unique[qobj$significant == TRUE]
sqtl_replicated$qvalue <- qobj$qvalues[qobj$significant == TRUE]

cat("\n  ✓ Replicated sQTLs (FDR <", FDR_THRESHOLD, "):",
    nrow(sqtl_replicated), "\n")
cat("  ✓ Unique genes:", uniqueN(sqtl_replicated$gene_id), "\n")
cat("  ✓ Replication rate:",
    round(100 * nrow(sqtl_replicated) / nrow(sqtl_trqtl_unique), 1), "%\n")

# Save results
fwrite(sqtl_replicated,
       sprintf("results/05_GTEx_sQTL_replicated_in_nanopore_FDR%.0f.txt",
               FDR_THRESHOLD * 100),
       sep = "\t")

cat("\n")

################################################################################
# PART 6: STRUCTURAL CATEGORY ANALYSIS (SQANTI3)
################################################################################

cat("PART 6: Comparing structural categories (SQANTI3)\n")
cat("----------------------------------------------------------------------\n")

# Load SQANTI3 classifications
# Expected location: merged_transcripts_flair_collapse_classification.txt.gz
sqanti <- fread("merged_transcripts_flair_collapse_classification.txt.gz")

cat("  ✓ Loaded SQANTI3 classifications:", nrow(sqanti), "transcripts\n")

# Load all significant nanopore trQTLs (for comparison)
# Expected location: LCLs_nanopore_FDR5_0PCs.significant.txt
nanopore_trqtl_sig <- fread("LCLs_nanopore_FDR5_0PCs.significant.txt")

# Clean column names
setnames(nanopore_trqtl_sig,
         old = c("V1", "V6"),
         new = c("gene_id", "transcript_id"),
         skip_absent = TRUE)

nanopore_trqtl_sig$gene_id_clean <- sub("\\..*", "", nanopore_trqtl_sig$gene_id)

cat("  ✓ Significant nanopore trQTLs:", nrow(nanopore_trqtl_sig), "\n")
cat("  ✓ Unique genes:", uniqueN(nanopore_trqtl_sig$gene_id_clean), "\n")

# Tag trQTLs as overlapping with sQTLs or not
genes_sqtl_replicated <- unique(sqtl_replicated$gene_id)

nanopore_trqtl_tagged <- nanopore_trqtl_sig %>%
  mutate(overlap_sqtl = gene_id_clean %in% genes_sqtl_replicated)

cat("  ✓ trQTLs overlapping with replicated sQTLs:",
    sum(nanopore_trqtl_tagged$overlap_sqtl), "\n")
cat("  ✓ trQTLs not overlapping (trQTL-only):",
    sum(!nanopore_trqtl_tagged$overlap_sqtl), "\n")

# Merge with SQANTI3 classifications
nanopore_with_sqanti <- nanopore_trqtl_tagged %>%
  left_join(
    sqanti %>% select(isoform, structural_category),
    by = c("transcript_id" = "isoform")
  ) %>%
  mutate(
    structural_category = ifelse(
      is.na(structural_category),
      "unclassified",
      trimws(structural_category)
    )
  )

# Summarize categories
cat("\n  Structural categories:\n")
print(table(nanopore_with_sqanti$structural_category))

# Test for enrichment of novel categories in overlap vs trQTL-only
nanopore_with_sqanti <- nanopore_with_sqanti %>%
  mutate(
    is_novel = structural_category %in% c("novel_in_catalog", "novel_not_in_catalog")
  )

contingency_table <- table(
  nanopore_with_sqanti$overlap_sqtl,
  nanopore_with_sqanti$is_novel
)

cat("\n  Contingency table (overlap × novel):\n")
print(contingency_table)

fisher_result <- fisher.test(contingency_table)
cat("\n  Fisher's exact test:\n")
cat("    Odds ratio:", round(fisher_result$estimate, 3), "\n")
cat("    P-value:", format.pval(fisher_result$p.value, digits = 3), "\n")

# Save analysis data
fwrite(nanopore_with_sqanti,
       "results/06_trQTL_with_SQANTI3_categories.txt",
       sep = "\t")

cat("\n")

################################################################################
# PART 7: GENERATE PLOTS
################################################################################

cat("PART 7: Generating plots\n")
cat("----------------------------------------------------------------------\n")

# Colorblind-friendly palette
color_fsm <- "#0072B2"    # Blue
color_novel <- "#E69F00"  # Orange
color_other <- "grey70"

# Prepare data for plotting
plot_data <- nanopore_with_sqanti %>%
  mutate(
    group = ifelse(overlap_sqtl, "sQTL ∩ trQTL", "trQTL only"),
    struct_simple = case_when(
      structural_category %in% c("novel_in_catalog", "novel_not_in_catalog") ~ "NIC/NNC",
      structural_category %in% c("full-splice_match", "incomplete-splice_match") ~ "FSM/ISM",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(group, struct_simple) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))

# Stacked bar plot
p1 <- ggplot(plot_data, aes(x = group, y = prop, fill = struct_simple)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = c(
      "FSM/ISM" = color_fsm,
      "NIC/NNC" = color_novel,
      "Other" = color_other
    ),
    name = "Category"
  ) +
  labs(
    x = "",
    y = "Fraction of transcripts",
    title = "SQANTI3 Structural Categories: sQTL-trQTL Overlap"
  ) +
  annotate(
    "text",
    x = 1.5, y = 1.05,
    label = sprintf("NIC/NNC depletion in overlap\nOR = %.2f, p = %.3f",
                   fisher_result$estimate,
                   fisher_result$p.value),
    size = 4
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave("plots/01_SQANTI3_sQTL_trQTL_structural_categories.pdf",
       plot = p1,
       width = 6,
       height = 5)

# NIC/NNC proportion comparison
nic_prop <- nanopore_with_sqanti %>%
  mutate(group = ifelse(overlap_sqtl, "sQTL ∩ trQTL", "trQTL only")) %>%
  group_by(group) %>%
  summarise(prop_novel = mean(is_novel))

p2 <- ggplot(nic_prop, aes(x = group, y = prop_novel, fill = group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(
    "sQTL ∩ trQTL" = color_fsm,
    "trQTL only" = "#009E73"
  )) +
  labs(
    x = "",
    y = "Fraction of NIC/NNC transcripts",
    title = "Novel Splice Junctions: Depletion in sQTL-trQTL Overlap"
  ) +
  annotate(
    "text",
    x = 1.5, y = max(nic_prop$prop_novel) * 1.1,
    label = sprintf("OR = %.2f, p = %.3f",
                   fisher_result$estimate,
                   fisher_result$p.value),
    size = 4
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave("plots/02_NIC_NNC_depletion_in_overlap.pdf",
       plot = p2,
       width = 6,
       height = 5)

cat("  ✓ Plots saved to plots/\n\n")

################################################################################
# PART 8: SUMMARY STATISTICS
################################################################################

cat("PART 8: Writing summary statistics\n")
cat("----------------------------------------------------------------------\n")

summary_stats <- data.table(
  Metric = c(
    "GTEx sQTL significant pairs (discovery)",
    "GTEx sQTL genes tested in nanopore",
    "Nanopore gene-SNP pairs tested",
    "Nanopore sQTLs replicated (FDR < 5%)",
    "Replication rate (%)",
    "trQTLs overlapping with replicated sQTLs",
    "trQTLs not overlapping (trQTL-only)",
    "Novel transcripts in overlap (NIC/NNC)",
    "Novel transcripts in trQTL-only (NIC/NNC)",
    "Fisher's OR (novel enrichment)",
    "Fisher's p-value"
  ),
  Value = c(
    nrow(gtex_sqtl_pairs),
    length(genes_common),
    nrow(sqtl_trqtl_unique),
    nrow(sqtl_replicated),
    round(100 * nrow(sqtl_replicated) / nrow(sqtl_trqtl_unique), 1),
    sum(nanopore_trqtl_tagged$overlap_sqtl),
    sum(!nanopore_trqtl_tagged$overlap_sqtl),
    sum(nanopore_with_sqanti$overlap_sqtl & nanopore_with_sqanti$is_novel),
    sum(!nanopore_with_sqanti$overlap_sqtl & nanopore_with_sqanti$is_novel),
    round(fisher_result$estimate, 3),
    format.pval(fisher_result$p.value, digits = 3)
  )
)

fwrite(summary_stats, "results/99_sQTL_replication_summary.txt", sep = "\t")

print(summary_stats)

cat("\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("======================================================================\n")
cat("sQTL Replication Analysis Complete!\n")
cat("======================================================================\n\n")

cat("Output files:\n")
cat("  results/01_sQTL_v37_v38_coordinates_with_rsIDs.txt\n")
cat("  results/02_GTEx_sQTL_complete_with_coordinates.txt\n")
cat("  results/03_GTEx_sQTL_testable_in_nanopore.txt\n")
cat("  results/04_GTEx_gene_SNP_pairs_for_filtering.txt\n")
cat("  results/05_GTEx_sQTL_replicated_in_nanopore_FDR5.txt\n")
cat("  results/06_trQTL_with_SQANTI3_categories.txt\n")
cat("  results/99_sQTL_replication_summary.txt\n")
cat("  plots/01_SQANTI3_sQTL_trQTL_structural_categories.pdf\n")
cat("  plots/02_NIC_NNC_depletion_in_overlap.pdf\n\n")

cat("Key findings:\n")
cat("  • sQTL replication rate:",
    round(100 * nrow(sqtl_replicated) / nrow(sqtl_trqtl_unique), 1), "%\n")
cat("  • π₁ (true signals):", round(1 - qobj$pi0, 3), "\n")
cat("  • Novel isoforms (NIC/NNC) depleted in overlap\n")
cat("    OR =", round(fisher_result$estimate, 2),
    "p =", format.pval(fisher_result$p.value, digits = 2), "\n")

cat("\n✓ Analysis complete!\n\n")
