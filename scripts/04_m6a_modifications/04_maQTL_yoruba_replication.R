#!/usr/bin/env Rscript
# =============================================================================
# maQTL Pipeline: Yoruba m6A Replication
# =============================================================================
# Description:
#   Replicates maQTL findings in an independent Yoruba LCL m6A-QTL dataset
#   (m6A-seq, antibody-based enrichment, 1000 Genomes LCLs).
##
#   Steps:
#   1. Genomic overlap — all DRS modifications vs Yoruba m6A-seq peaks,
#        by (a) position only and (b) position + same gene symbol.
#
#   2. Overlap restricted to m6A DRS modifications tested for maQTL.
#
#   3. Significant maQTL hits vs Yoruba peaks + Yoruba m6A-QTL.
#
#   4. UpSet plot — set membership across all / m6aQTL-tested / m6aQTL-significant sets
#        relative to Yoruba peak overlap.
#
# Input:
#   m6AQTL.m6APeak_logOR_GC.IP.adjusted_qqnorm.15PCs.fastQTL.nominals.rds
#                                          (Yoruba m6A-QTL nominals, all peaks)
#   m6APeak_logOR_GC.IP.adjusted_qqnorm.fastQTL.txt.gz
#                                          (Yoruba m6A-seq peak coordinates)
#   filtered_m6a_modifications.txt            (Script 01 — filtered site table with
#                                              genomic coords: chr, start, end,
#                                              transcript, position, gene_id, gene_symbol)
#   LCLs_m6A_FDR10_3PCs.significant.txt      (Script 02 — significant maQTL hits)
#
# Output:
#   overlap_summary_stats.txt
#   nano_overlap_peaks_by_position.txt
#   nano_overlap_peaks_by_position_AND_gene.txt
#   maQTL_hits_yoruba_peak_status.txt
#   CALR_replication_yoruba.txt
#   upset_modifications_overlap_maQTL_tested_sig.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(UpSetR)
})

# =============================================================================
# CONFIG
# =============================================================================

YORUBA_RDS     <- "m6AQTL.m6APeak_logOR_GC.IP.adjusted_qqnorm.15PCs.fastQTL.nominals.rds"
PEAK_FILE      <- "m6APeak_logOR_GC.IP.adjusted_qqnorm.fastQTL.txt.gz"
FILTERED_MODS  <- "filtered_m6a_modifications.txt"  # Script 01 output
MAQTL_SIG      <- "LCLs_m6A_FDR10_3PCs.significant.txt"  # Script 02 output

# CALR lead SNP from the maQTL analysis
CALR_LEAD_SNP  <- "rs2974751"

# =============================================================================
# HELPERS
# =============================================================================

log_step <- function(...) message("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

pct <- function(n, d) paste0(round(100 * n / d, 1), "%")

# Build GRanges from filtered_m6a_modifications.txt rows
# Required columns: chr, modification_start (or start), modification_end (or end),
#                   transcript (transcript_id), gene_id, gene_symbol
mods_to_gr <- function(dt) {
  # Handle both 'start'/'end' and 'modification_start'/'modification_end' column names
  s <- if ("modification_start" %in% names(dt)) dt$modification_start else dt$start
  e <- if ("modification_end"   %in% names(dt)) dt$modification_end   else dt$end
  chr_col <- if ("chr" %in% names(dt)) dt$chr else dt[[1]]
  GRanges(
    seqnames    = chr_col,
    ranges      = IRanges(start = as.integer(s), end = as.integer(e)),
    mod_id      = dt$modification_id,
    transcript  = dt$transcript,
    gene_ensg   = dt$gene_id,
    gene_symbol = dt$gene_symbol
  )
}

# =============================================================================
# LOAD YORUBA PEAKS
# =============================================================================

log_step("Loading Yoruba m6A-seq peaks...")

if (!file.exists(PEAK_FILE)) stop("Peak file not found: ", PEAK_FILE)

# PEAK column format: chr12:113000000-113010000_OAS1_+
peak_raw <- fread(PEAK_FILE)
setDT(peak_raw)
peak_raw[, peak_chr    := sub(":.*",                  "", PEAK)]
peak_raw[, peak_start  := as.integer(sub(".*:(\\d+)-.*",   "\\1", PEAK))]
peak_raw[, peak_end    := as.integer(sub(".*-(\\d+)_.*",   "\\1", PEAK))]
peak_raw[, peak_gene   := sub(".*_(.*)_[+-]$",        "\\1", PEAK)]
peak_raw[, peak_strand := sub(".*_([+-])$",            "\\1", PEAK)]
message("  Yoruba m6A-seq peaks: ", nrow(peak_raw),
        " | unique genes: ", uniqueN(peak_raw$peak_gene))

gr_peak <- GRanges(
  seqnames  = peak_raw$peak_chr,
  ranges    = IRanges(start = peak_raw$peak_start, end = peak_raw$peak_end),
  PEAK      = peak_raw$PEAK,
  gene_peak = peak_raw$peak_gene
)



# Helper: run both overlap flavours and return summary
run_overlap <- function(gr_mods, label) {
  hits   <- findOverlaps(gr_mods, gr_peak, ignore.strand = TRUE)
  ol <- data.table(
    mod_id      = mcols(gr_mods)$mod_id[queryHits(hits)],
    gene_symbol = mcols(gr_mods)$gene_symbol[queryHits(hits)],
    gene_ensg   = mcols(gr_mods)$gene_ensg[queryHits(hits)],
    transcript  = mcols(gr_mods)$transcript[queryHits(hits)],
    chr         = as.character(seqnames(gr_mods))[queryHits(hits)],
    mod_start   = start(gr_mods)[queryHits(hits)],
    peak        = mcols(gr_peak)$PEAK[subjectHits(hits)],
    gene_peak   = mcols(gr_peak)$gene_peak[subjectHits(hits)]
  )

  # Position only
  n_pos  <- uniqueN(ol$mod_id)
  # Position + same gene
  ol_g   <- ol[!is.na(gene_symbol) & !is.na(gene_peak) &
                 toupper(trimws(gene_symbol)) == toupper(trimws(gene_peak))]
  n_gene <- uniqueN(ol_g$mod_id)
  n_tot  <- length(unique(mcols(gr_mods)$mod_id))

  message("  [", label, "] n=", n_tot)
  message("    Position only      : ", n_pos,  " (", pct(n_pos,  n_tot), ")")
  message("    Position + gene    : ", n_gene, " (", pct(n_gene, n_tot), ")")

  list(ol_pos = ol, ol_gene = ol_g, n_tot = n_tot,
       n_pos = n_pos, n_gene = n_gene)
}

# =============================================================================
# STEP 1: All DRS modifications vs Yoruba peaks
# =============================================================================

log_step("Step 1+2: Loading filtered modifications and overlapping with Yoruba peaks...")

# Load filtered modifications from Script 01
# Columns: transcript, position, gene_id, chr, start, end, strand, source, [gene_symbol]
filt <- fread(FILTERED_MODS, header = TRUE, sep = "\t")
message("  Filtered modifications loaded: ", nrow(filt))

# Build modification_id and gene_symbol if not already present
if (!"modification_id" %in% names(filt))
  filt[, modification_id := paste0(transcript, "_", position)]
if (!"gene_symbol" %in% names(filt))
  filt[, gene_symbol := gene_id]   # replace with actual symbol column if available

# ALL mods = all rows in filtered_m6a_modifications.txt
gr_all  <- mods_to_gr(filt)
res_all <- run_overlap(gr_all, "ALL mods")

fwrite(res_all$ol_pos,  "nano_overlap_peaks_by_position.txt",          sep = "\t")
fwrite(res_all$ol_gene, "nano_overlap_peaks_by_position_AND_gene.txt", sep = "\t")
message("  Written: nano_overlap_peaks_by_position*.txt")

# maQTL-tested mods = same set (≥50% filter already applied in Script 01)
# Run separately so stats are reported independently in the summary table
gr_mod  <- gr_all
res_mod <- run_overlap(gr_mod, "maQTL-tested mods")

# =============================================================================
# STEP 3: Significant maQTL hits — peak overlap + Yoruba QTL replication
# =============================================================================

log_step("Step 3: Significant maQTL hits — Yoruba peak and QTL overlap...")

maqtl_sig  <- fread(MAQTL_SIG, header = FALSE)
# V1=transcript_id  V6=gene_id  V10=lead_SNP  (QTLtools permutation output)
maqtl_dt   <- data.table(
  transcript_id = maqtl_sig$V1,
  gene_id       = sub("\\..*", "", maqtl_sig$V6),
  gene_symbol   = maqtl_sig$V7,   # adjust column if needed
  lead_snp      = maqtl_sig$V10,
  beta          = maqtl_sig$V14,
  pvalue        = maqtl_sig$V18
)
message("  Significant maQTL hits: ", nrow(maqtl_dt))

# Subset filtered mods to significant maQTL hits using MAQTL_SIG
# V1 in significant file = pid (transcript_id + "_" + position + "_" + motif)
sig_dt  <- fread(MAQTL_SIG, header = FALSE)
sig_ids <- unique(sub("_[ACGT]{5}$", "", sig_dt$V1))  # strip motif suffix if present
filt[, mod_base_id := paste0(transcript, "_", position)]
nano_maqtl <- filt[mod_base_id %in% sig_ids]
message("  Significant maQTL modifications: ", nrow(nano_maqtl))
gr_maqtl   <- mods_to_gr(nano_maqtl)
res_maqtl  <- run_overlap(gr_maqtl, "maQTL sig mods")

mods_in_peak <- unique(res_maqtl$ol_gene$mod_id)
mods_maqtl   <- unique(mcols(gr_maqtl)$mod_id)

# Flag each maQTL hit by whether its modification overlaps a Yoruba peak
maqtl_peak_status <- data.table(
  mod_id      = mods_maqtl,
  in_peak     = mods_maqtl %in% mods_in_peak
)
maqtl_peak_status <- merge(
  maqtl_peak_status,
  nano_maqtl[, .(modification_id, gene_symbol, transcript_id = transcript)],
  by.x = "mod_id", by.y = "modification_id",
  all.x = TRUE
)

message("\n  maQTL modifications overlapping a same-gene Yoruba peak: ",
        sum(maqtl_peak_status$in_peak), " / ", nrow(maqtl_peak_status))

not_in_peak <- maqtl_peak_status[in_peak == FALSE]
message("  maQTL modifications with NO Yoruba peak: ", nrow(not_in_peak))
if (nrow(not_in_peak) > 0) {
  message("  Gene(s) without peak:")
  print(not_in_peak[, .(mod_id, gene_symbol, transcript_id)])
}

fwrite(maqtl_peak_status, "maQTL_hits_yoruba_peak_status.txt", sep = "\t")
message("  Written: maQTL_hits_yoruba_peak_status.txt")

# =============================================================================
# STEP 4: Yoruba m6A-QTL replication — lead SNP + direction of effect
# =============================================================================

log_step("Step 4: Yoruba m6A-QTL replication — lead SNP and effect direction...")

if (!file.exists(YORUBA_RDS)) {
  message("  WARNING: Yoruba QTL RDS not found: ", YORUBA_RDS, " — skipping")
} else {

  m6a_all <- readRDS(YORUBA_RDS)
  setDT(m6a_all)
  m6a_sig <- m6a_all[qvalue < 0.05]
  message("  Yoruba FDR5 m6A-QTL sites: ", nrow(m6a_sig),
          " | unique genes: ", uniqueN(m6a_sig$GENE))

  # Match each maQTL lead SNP against Yoruba significant sites
  # Join on gene symbol + SNP ID
  qtl_match <- merge(
    maqtl_dt[, .(gene_symbol, lead_snp, beta, pvalue)],
    m6a_sig[,  .(GENE, SNPID, beta_yoruba = beta, pval_yoruba = pvalue, qval_yoruba = qvalue)],
    by.x = c("gene_symbol", "lead_snp"),
    by.y = c("GENE",        "SNPID"),
    all.x = FALSE
  )

  message("\n  maQTL lead SNPs replicated in Yoruba m6A-QTL: ",
          nrow(qtl_match), " / ", nrow(maqtl_dt))

  if (nrow(qtl_match) > 0) {
    # Direction of effect: same sign = concordant
    qtl_match[, direction_concordant := sign(beta) == sign(beta_yoruba)]
    message("  Direction concordant: ", sum(qtl_match$direction_concordant),
            " / ", nrow(qtl_match))
    print(qtl_match[, .(gene_symbol, lead_snp, beta, beta_yoruba,
                         pvalue, pval_yoruba, direction_concordant)])
  }

  # CALR specific check
  calr_yoruba <- m6a_sig[GENE == "CALR" & SNPID == CALR_LEAD_SNP]
  calr_maqtl  <- maqtl_dt[gene_symbol == "CALR"]

  message("\n  CALR replication check (lead SNP: ", CALR_LEAD_SNP, "):")
  if (nrow(calr_yoruba) > 0 && nrow(calr_maqtl) > 0) {
    calr_out <- data.table(
      gene             = "CALR",
      lead_snp         = CALR_LEAD_SNP,
      beta_maQTL       = calr_maqtl$beta[1],
      pval_maQTL       = calr_maqtl$pvalue[1],
      beta_Yoruba      = calr_yoruba$beta[1],
      pval_Yoruba      = calr_yoruba$pvalue[1],
      qval_Yoruba      = calr_yoruba$qvalue[1],
      direction_concordant = sign(calr_maqtl$beta[1]) == sign(calr_yoruba$beta[1])
    )
    print(calr_out)
    fwrite(calr_out, "CALR_replication_yoruba.txt", sep = "\t")
    message("  Written: CALR_replication_yoruba.txt")
  } else {
    message("  CALR or lead SNP ", CALR_LEAD_SNP, " not found in Yoruba FDR5 sites")
    message("  Available CALR entries in Yoruba (all q-values):")
    print(m6a_all[GENE == "CALR" & SNPID == CALR_LEAD_SNP,
                  .(GENE, SNPID, beta, pvalue, qvalue)])
  }
}

# =============================================================================
# STEP 5: Summary Statistics Table (paper numbers)
# =============================================================================

log_step("Step 5: Summary statistics...")

summary_stats <- data.table(
  Set                   = c("All DRS mods",
                             "All DRS mods",
                             "maQTL-tested mods (n=2,761)",
                             "Significant maQTL hits"),
  Overlap_type          = c("Position only (relaxed)",
                             "Position + same gene",
                             "Position + same gene",
                             "Position + same gene"),
  N_mods                = c(res_all$n_tot,   res_all$n_tot,
                             res_mod$n_tot, length(mods_maqtl)),
  N_overlapping         = c(res_all$n_pos,   res_all$n_gene,
                             res_mod$n_gene, length(mods_in_peak)),
  Pct_overlapping       = c(pct(res_all$n_pos,   res_all$n_tot),
                             pct(res_all$n_gene,  res_all$n_tot),
                             pct(res_mod$n_gene, res_mod$n_tot),
                             pct(length(mods_in_peak), length(mods_maqtl)))
)

message("\n  ── Overlap summary (paper statistics) ──")
print(summary_stats)
fwrite(summary_stats, "overlap_summary_stats.txt", sep = "\t")
message("  Written: overlap_summary_stats.txt")

# =============================================================================
# STEP 6: UpSet Plot — all / maQTL-tested / maQTL-significant sets vs Yoruba peak overlap
# =============================================================================

log_step("Step 6: UpSet plot...")

mods_universe   <- unique(filt$modification_id)
mods_tested_ids <- mods_universe                      # same set as ALL
mods_maqtl_ids  <- unique(nano_maqtl$modification_id)
mods_pos_peak  <- unique(res_all$ol_pos$mod_id)   # position-only overlap
mods_gene_peak <- unique(res_all$ol_gene$mod_id)  # position + gene overlap

upset_df <- data.frame(
  Yoruba_peak_sameGene = as.integer(mods_universe %in% mods_gene_peak),
  maQTL_tested         = as.integer(mods_universe %in% mods_tested_ids),
  maQTL_significant    = as.integer(mods_universe %in% mods_maqtl_ids)
)

message("  Set sizes:")
message("    Yoruba peak (same gene): ", sum(upset_df$Yoruba_peak_sameGene))
message("    maQTL tested           : ", sum(upset_df$maQTL_tested))
message("    maQTL significant      : ", sum(upset_df$maQTL_significant))

pdf("upset_modifications_overlap_maQTL_tested_sig.pdf", width = 10, height = 6)
UpSetR::upset(
  upset_df,
  sets            = c("Yoruba_peak_sameGene", "maQTL_tested", "maQTL_significant"),
  keep.order      = TRUE,
  order.by        = "freq",
  mainbar.y.label = "Intersection size",
  sets.x.label    = "Set size",
  text.scale      = 1.4
)
dev.off()
message("  Written: upset_modifications_overlap_maQTL_tested_sig.pdf")

message("\nScript 03 complete.")
message("Next: Rscript 04_maQTL_visualization.R")
