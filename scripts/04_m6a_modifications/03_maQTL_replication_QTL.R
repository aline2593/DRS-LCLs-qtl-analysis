#!/usr/bin/env Rscript
# =============================================================================
# maQTL Pipeline: QTL Replication (eQTL / sQTL)
# =============================================================================
# Description:
#   Tests whether maQTL lead variants replicate in two independent short-read
#   QTL datasets using the nominal pass results:
#
#   1. Illumina eQTL — gene-level expression QTLs (317 samples LCLs)
#   2. GTEx sQTL     — splicing QTLs in EBV-transformed lymphocytes GTEx (LCLs)
#
#   For each dataset:
#     - Matches maQTL nominal results to external QTL significant hits
#       on (gene_id, SNP_id)
#     - Applies q-value FDR correction to the maQTL nominal p-values
#       within the overlapping set
#     - Writes significant overlapping hits
#     - Prints a cross-dataset summary table
#
# Input:
#   Nominal_gene_ID.txt                                         (Script 01)
#   LCLs_m6A_FDR10_3PCs.significant.txt                         (Script 01)
#   LCL_RNA.chunkALL.significant_permutation.txt.gz             (Illumina eQTL)
#   Cells_EBV-transformed_lymphocytes.v8.sgenes_mod.txt         (GTEx sQTL)
#
# Output:
#   Illumina_eQTL_extracted.txt
#   eQTL_recap_modnominal_nanopore.txt
#   FDR5_eQTLs_gene_recap_modifications.txt
#   sQTL_recap_modnominal_nanopore.txt
#   FDR5_sQTLs_gene_recap_modifications.txt
#   QTL_replication_summary.txt
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(qvalue)
})

# =============================================================================
# CONFIG
# =============================================================================

NOMINAL_FILE  <- "Nominal_gene_ID.txt"
MAQTL_SIG     <- "LCLs_m6A_FDR10_3PCs.significant.txt"

EQTL_FILE  <- "LCL_RNA.chunkALL.significant_permutation.txt.gz"
SQTL_FILE  <- "Cells_EBV-transformed_lymphocytes.v8.sgenes_mod.txt"

# =============================================================================
# HELPERS
# =============================================================================

log_step <- function(...) message("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

strip_version <- function(x) sub("\\..*", "", x)

# Apply q-value FDR; fall back to BH when too few p-values
apply_qvalue_fdr <- function(df, pval_col, fdr_threshold = 0.05) {
  pvals <- df[[pval_col]]
  if (sum(!is.na(pvals)) < 10) {
    message("  Too few p-values for qvalue — using BH correction")
    df$q_adj <- p.adjust(pvals, method = "BH")
    return(df[!is.na(df$q_adj) & df$q_adj <= fdr_threshold, , drop = FALSE])
  }
  q          <- qvalue(pvals, fdr.level = fdr_threshold)
  df$q_value <- q$qvalues
  summary(q)
  df[q$significant, , drop = FALSE]
}

# =============================================================================
# LOAD NOMINAL RESULTS
# =============================================================================

log_step("Loading nominal maQTL results...")

# Columns: V1=mod_id  V2=SNP  V3=nominal_pvalue  V4=gene_id
Nominal <- fread(NOMINAL_FILE, header = FALSE, sep = " ")
Nominal[, V4 := strip_version(V4)]
message("  Nominal: ", nrow(Nominal), " gene–SNP pairs")

# =============================================================================
# STEP 1: eQTL Replication — Illumina gene expression
# =============================================================================

log_step("Step 1: eQTL replication (Illumina gene expression)...")

eqtl_raw <- fread(EQTL_FILE, header = FALSE, sep = " ")

# V1=gene_id  V8=SNP_id  V16=nominal_pvalue
eqtl <- eqtl_raw[, .(
  gene_id      = strip_version(V1),
  SNP_id       = V8,
  eQTL_nominal = V16
)]
fwrite(eqtl, "Illumina_eQTL_extracted.txt", sep = "\t")
message("  Illumina eQTL hits: ", nrow(eqtl))

# Match on (gene_id, SNP_id)
tot_eqtl <- merge(
  Nominal, eqtl,
  by.x = c("V4", "V2"), by.y = c("gene_id", "SNP_id")
)
message("  Overlapping gene-SNP pairs: ", nrow(tot_eqtl),
        " | unique genes: ", uniqueN(tot_eqtl$V4))

# One row per gene: keep strongest maQTL nominal p-value
tot_eqtl_uniq <- tot_eqtl[order(V3)][!duplicated(V4)]
setnames(tot_eqtl_uniq,
  c("V4", "V2", "V1", "V3", "eQTL_nominal"),
  c("Gene_id", "SNP_ID", "Gene_mod_ID", "nominal_pvalue_mod", "eQTL_nominal")
)
fwrite(tot_eqtl_uniq, "eQTL_recap_modnominal_nanopore.txt", sep = "\t")

sig_eqtl <- apply_qvalue_fdr(as.data.frame(tot_eqtl_uniq), "nominal_pvalue_mod")
message("  Significant at FDR 5%: ", nrow(sig_eqtl))
fwrite(as.data.table(sig_eqtl), "FDR5_eQTLs_gene_recap_modifications.txt", sep = "\t")

if (nrow(sig_eqtl) > 0) {
  message("  Significant hits:")
  print(sig_eqtl[, c("Gene_id", "SNP_ID", "Gene_mod_ID",
                      "nominal_pvalue_mod", "eQTL_nominal")])
}

# =============================================================================
# STEP 2: sQTL Replication — GTEx EBV-transformed lymphocytes
# =============================================================================

log_step("Step 2: sQTL replication (GTEx EBV-transformed lymphocytes)...")

sqtl_raw <- fread(SQTL_FILE, header = TRUE, sep = "\t")

# Keep: gene_id, rsID, nominal p-value, alleles, TSS distance
sqtl <- sqtl_raw[, .(
  gene_id      = strip_version(gene_id),
  SNP_id       = rs_id_dbSNP151_GRCh38p7,
  sQTL_nominal = pval_nominal,
  ref          = ref,
  alt          = alt,
  tss_distance = tss_distance
)]
message("  GTEx sQTL hits: ", nrow(sqtl))

# Work from a clean copy to avoid overwriting V4 strip already done above
Nominal_s <- copy(Nominal)
Nominal_s[, V4 := strip_version(V4)]

tot_sqtl <- merge(
  Nominal_s, sqtl,
  by.x = c("V4", "V2"), by.y = c("gene_id", "SNP_id")
)
message("  Overlapping gene-SNP pairs: ", nrow(tot_sqtl),
        " | unique genes: ", uniqueN(tot_sqtl$V4))

tot_sqtl_uniq <- tot_sqtl[order(V3)][!duplicated(V4)]
setnames(tot_sqtl_uniq,
  c("V4", "V2", "V1", "V3"),
  c("Gene_id", "SNP_ID", "Gene_mod_ID", "nominal_pvalue_mod")
)
fwrite(tot_sqtl_uniq, "sQTL_recap_modnominal_nanopore.txt", sep = "\t")

sig_sqtl <- apply_qvalue_fdr(as.data.frame(tot_sqtl_uniq), "nominal_pvalue_mod")
message("  Significant at FDR 5%: ", nrow(sig_sqtl))
fwrite(as.data.table(sig_sqtl), "FDR5_sQTLs_gene_recap_modifications.txt", sep = "\t")

if (nrow(sig_sqtl) > 0) {
  message("  Significant hits:")
  print(sig_sqtl[, c("Gene_id", "SNP_ID", "Gene_mod_ID",
                      "nominal_pvalue_mod", "sQTL_nominal",
                      "ref", "alt", "tss_distance")])
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

log_step("Summary: cross-dataset QTL replication...")

summary_tbl <- data.table(
  Dataset           = c("Illumina eQTL", "GTEx sQTL"),
  Overlapping_pairs = c(nrow(tot_eqtl_uniq), nrow(tot_sqtl_uniq)),
  Unique_genes      = c(uniqueN(tot_eqtl_uniq$Gene_id),
                        uniqueN(tot_sqtl_uniq$Gene_id)),
  Sig_FDR5          = c(nrow(sig_eqtl), nrow(sig_sqtl))
)

print(summary_tbl)
fwrite(summary_tbl, "QTL_replication_summary.txt", sep = "\t")
message("  Summary → QTL_replication_summary.txt")

message("\nScript 02 complete.")
message("Next: Rscript 03_maQTL_yoruba_replication.R")
