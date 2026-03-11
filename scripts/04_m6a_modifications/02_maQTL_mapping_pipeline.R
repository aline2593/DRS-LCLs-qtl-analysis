#!/usr/bin/env Rscript
# =============================================================================
# maQTL Pipeline: Full QTL Mapping
# =============================================================================
# Description:
#   A — BED Preparation
#         Load Nanopore_ModificationsRatio.txt (per-sample modification ratios),
#         restrict to sites that passed filtering in Script 01
#         (filtered_m6a_modifications.txt), attach gene-level coordinates
#         from Ref_file_readytouse.txt, rename GM→NA sample IDs, sort,
#         bgzip, tabix index.
#
#   B — QTLtools PCA + Covariate Files
#         Run QTLtools PCA on phenotypes, merge with fixed covariates
#         (SEX, GEN, GPC1-3), create merged.pcs.{0,3,6,...,60}.
#
#   C — Permutation Pass
#         Submit QTLtools cis --permute jobs for all PC levels (0–60, step 3)
#         across 20 chunks via LSF. After jobs complete, merge chunks.
#
#   D — FDR Calibration & PC Selection
#         Compute q-values per PC level, plot FDR and pi1 curves,
#         extract significant hits at the selected PC level.
#
#   E — Nominal Pass
#         Prepare BED with composite pid_gid ID, submit 20-chunk nominal
#         jobs, merge chunks, extract key columns, join with gene IDs.
#
# Run modes:
#   Rscript 02_maQTL_mapping_pipeline.R                  # Parts A–B + submit C
#   Rscript 02_maQTL_mapping_pipeline.R --merge-perm     # merge Part C chunks
#   Rscript 02_maQTL_mapping_pipeline.R --fdr            # Part D
#   Rscript 02_maQTL_mapping_pipeline.R --nominal        # submit Part E
#   Rscript 02_maQTL_mapping_pipeline.R --merge-nominal  # merge Part E chunks
#
# Input:
#   Nanopore_ModificationsRatio.txt      (Script 00 — per-sample modification ratios)
#   filtered_m6a_modifications.txt       (Script 01 — filtered site list)
#   Ref_file_readytouse.txt              (Script 01 — gene annotation coordinates)
#   merged.pcs.0   (fixed covariates: SEX, GEN, GPC1-3 — pre-existing)
#
# Output:
#   Nanopore_ModificationsRatio_NEW.bed.gz (.tbi)
#   genes.50percent.pca
#   merged.pcs  /  merged.pcs.{0,3,...,60}
#   m6A_perm1000_cov{N}_merged.txt.gz
#   pval_dists.pdf  /  no_sig_by_pcs_FDRs.pdf
#   LCLs_m6A_FDR{5,10}_{N}PCs.significant.txt
#   Nanopore_ModificationsRatio_nominal.bed.gz (.tbi)
#   Nominal_gene_ID.txt
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(qvalue)
})

args          <- commandArgs(trailingOnly = TRUE)
MERGE_PERM    <- "--merge-perm"    %in% args
FDR_MODE      <- "--fdr"           %in% args
NOMINAL_MODE  <- "--nominal"       %in% args
MERGE_NOMINAL <- "--merge-nominal" %in% args
RUN_PREP      <- !any(c(MERGE_PERM, FDR_MODE, NOMINAL_MODE, MERGE_NOMINAL))

# =============================================================================
# CONFIG — edit paths here
# =============================================================================

MOD_MATRIX    <- "Nanopore_ModificationsRatio.txt"     # Script 00 output
FILTERED_MODS <- "filtered_m6a_modifications.txt"    # Script 01 output
GENE_REF      <- "Ref_file_readytouse.txt"           # Script 01 output
OUTPUT_BED    <- "Nanopore_ModificationsRatio_NEW.bed"
OUTPUT_BED_GZ <- paste0(OUTPUT_BED, ".gz")
BED_NOMINAL   <- "Nanopore_ModificationsRatio_nominal.bed.gz"

VCF           <- "path/to/genotypes.bcf"
COV_BASE      <- "merged.pcs.0"       # pre-existing fixed covariates file
PC_PREFIX     <- "genes.50percent"    # QTLtools PCA output prefix
PC_MERGED     <- "merged.pcs"         # final merged covariate file

SELECTED_PC   <- 3                    # update after reviewing FDR plots
FDR_SCRIPT    <- "qtltools_runFRD_cis.R"
OUTPUT_PREFIX <- "LCLs_m6A"

N_CHUNKS      <- 20

CHR_ORDER     <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# =============================================================================
# HELPERS
# =============================================================================

log_step <- function(...) message("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

submit_lsf <- function(cmd, jobname, walltime = "20:10",
                       stdout = NULL, stderr = NULL) {
  out_flag <- if (!is.null(stdout)) paste("-o", stdout) else ""
  err_flag <- if (!is.null(stderr)) paste("-e", stderr) else ""
  bsub_cmd <- paste(
    "echo", shQuote(cmd), "|",
    "bsub",
    "-P", LSF_PROJECT,
    "-q", LSF_QUEUE,
    "-n 1",
    "-W", walltime,
    "-J", jobname,
    out_flag, err_flag,
    "-R", shQuote(paste0("rusage[mem=", LSF_MEM, "]"))
  )
  system(bsub_cmd)
}

# =============================================================================
# PART A: BED Preparation
# =============================================================================

if (RUN_PREP) {

  log_step("Part A: BED preparation...")

  # -- A1. Load the full per-sample modification ratio matrix (Script 00 output)
  # Columns: transcript_id, transcript_position, [site cols], id, [sample cols...]
  mod_raw <- fread(MOD_MATRIX, header = TRUE, sep = "\t")
  message("  Loaded modification matrix: ", nrow(mod_raw), " sites x ", ncol(mod_raw), " cols")

  # -- A2. Load filtered site list (Script 01 output) and restrict to passing sites
  # filtered_m6a_modifications.txt columns: transcript, position, gene_id, chr, ...
  filtered <- fread(FILTERED_MODS, header = TRUE, sep = "\t")
  message("  Filtered sites from Script 01: ", nrow(filtered))

  # Build site key matching the composite 'id' column in mod_raw
  filtered[, site_key := paste0(transcript, "_", position)]

  # Sample columns start at col 7 in mod_raw (one column per individual)
  sample_cols <- colnames(mod_raw)[7:ncol(mod_raw)]

  mod_filtered <- mod_raw[id %in% filtered$site_key]
  message("  Sites retained after filter join: ", nrow(mod_filtered))

  # -- A3. Load gene-level coordinates (Script 01 output)
  # Columns: chr, source, feature, start, end, score, strand, gene_id, transcript_id
  # Gene coordinates define the QTLtools cis-window (not the exact m6A position)
  ref <- fread(GENE_REF, header = TRUE, sep = "\t")
  gene_coords <- ref[feature == "transcript", .(
    chr    = chr[1],
    start  = min(as.integer(start)),
    end    = max(as.integer(end)),
    strand = strand[1]
  ), by = gene_id]

  # -- A4. Attach gene coordinates and gene_id to filtered sites
  mod_filtered[, pid := id]
  mod_filtered <- merge(
    mod_filtered,
    filtered[, .(site_key, gene_id)],
    by.x = "pid", by.y = "site_key",
    all.x = TRUE
  )
  mod_filtered <- merge(
    mod_filtered,
    gene_coords,
    by    = "gene_id",
    all.x = TRUE
  )
  message("  Sites with gene coordinates: ", nrow(mod_filtered))

  # -- A5. Build QTLtools BED: #Chr, start, end, pid, gid, strand, [samples...]
  mod_filtered[, `#Chr` := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  mod_filtered[, `#Chr` := factor(`#Chr`, levels = CHR_ORDER)]
  setorder(mod_filtered, `#Chr`, start, end)
  mod_filtered[, `#Chr` := as.character(`#Chr`)]

  col_order <- c("#Chr", "start", "end", "pid", "gene_id", "strand", sample_cols)
  col_order <- col_order[col_order %in% colnames(mod_filtered)]
  bed <- mod_filtered[, col_order, with = FALSE]
  setnames(bed, "gene_id", "gid")

  # Rename GM → NA in sample column names
  colnames(bed) <- gsub("GM", "NA", colnames(bed))

  message("  Final BED: ", nrow(bed), " sites x ", ncol(bed), " cols")

  write.table(bed, OUTPUT_BED, sep = "\t", row.names = FALSE, quote = FALSE)
  system(paste("bgzip -f", OUTPUT_BED))
  system(paste("tabix -p bed", OUTPUT_BED_GZ))
  message("  Output: ", OUTPUT_BED_GZ, " (tabix indexed)")

  # ============================================================================
  # PART B: QTLtools PCA + Covariate Files
  # ============================================================================

  log_step("Part B: QTLtools PCA...")

  pca_cmd <- paste(
    "QTLtools pca",
    "--bed", OUTPUT_BED_GZ,
    "--scale --center",
    "--out", PC_PREFIX
  )
  message("  Running: ", pca_cmd)
  system(paste("module load qtltools/1.3 &&", pca_cmd))

  # Merge PCA output with fixed covariates
  system(paste("cp", COV_BASE, "merged.pcs.0_filtered"))
  system("sed -i '7d' merged.pcs.0_filtered")
  system(paste("cp merged.pcs.0_filtered", PC_MERGED))
  system(paste("sed -i '6q'", PC_MERGED))
  system(paste("cat", paste0(PC_PREFIX, ".pca"), ">>", PC_MERGED))
  system(paste("sed -i '7d'", PC_MERGED))
  message("  Merged covariate file: ", PC_MERGED)

  # Create per-PC covariate files: merged.pcs.{0, 3, 6, ..., 60}
  log_step("Part B: Creating per-PC covariate files...")

  cov_full   <- fread(PC_MERGED)
  fixed_rows <- c("SEX", "GEN", "GPC1", "GPC2", "GPC3")
  pc_rows    <- grep("^genes.*PC[0-9]+$", cov_full$Sample, value = TRUE)

  for (i in seq(0, 60, 3)) {
    rows_to_keep <- if (i == 0) fixed_rows else c(fixed_rows, pc_rows[seq_len(i)])
    cov_out  <- cov_full[Sample %in% rows_to_keep]
    out_file <- paste0("merged.pcs.", i)
    fwrite(cov_out, out_file, sep = "\t", quote = FALSE)
    message("  Saved: ", out_file, " (", nrow(cov_out), " covariates)")
  }

  # ============================================================================
  # PART C: Submit Permutation Pass
  # ============================================================================

  log_step("Part C: Submitting permutation pass jobs...")
  system("module load qtltools/1.3")

  for (cov_num in seq(0, 60, 3)) {
    cov_file <- paste0("PC/merged.pcs.", cov_num)
    if (!file.exists(cov_file)) {
      message("  WARNING: missing ", cov_file, " — skipping")
      next
    }
    for (j in seq_len(N_CHUNKS)) {
      cmd <- paste(
        "QTLtools cis",
        "--vcf", VCF,
        "--bed", OUTPUT_BED_GZ,
        "--cov", cov_file,
        "--normal --permute 1000 --grp-best",
        "--chunk", j, N_CHUNKS,
        "--out", paste0("maQTL_perm1000_", cov_num, "_", j, "_", N_CHUNKS, ".txt")
      )
      submit_lsf(cmd,
                 jobname = paste0("perm_pc", cov_num, "_", j),
                 stdout  = paste0("maQTL_perm_norm_pc", cov_num, ".out"),
                 stderr  = paste0("maQTL_perm_norm_pc", cov_num, ".err"))
    }
    message("  Submitted ", N_CHUNKS, " jobs for PC=", cov_num)
  }

  message("\n  All permutation jobs submitted.")
  message("  Wait for completion, then run:")
  message("    Rscript 01_maQTL_mapping_pipeline.R --merge-perm")

} # end RUN_PREP

# =============================================================================
# PART C (cont.): Merge Permutation Chunks
# =============================================================================

if (MERGE_PERM) {

  log_step("Part C: Merging permutation chunks...")

  for (cov in seq(0, 60, 3)) {
    files <- list.files(pattern = paste0("^maQTL_perm1000_", cov, "_.*\\.txt$"))
    if (length(files) == 0) {
      message("  WARNING: no files for cov=", cov)
      next
    }
    out    <- paste0("m6a_perm1000_cov", cov, "_merged.txt.gz")
    merged <- rbindlist(lapply(files, fread, header = FALSE))
    fwrite(merged, out, sep = " ", col.names = FALSE, compress = "gzip")
    message("  cov=", cov, " → ", out, " (", nrow(merged), " rows)")
  }

  message("\n  Merging complete.")
  message("  Next: Rscript 01_maQTL_mapping_pipeline.R --fdr")

} # end MERGE_PERM

# =============================================================================
# PART D: FDR Calibration and PC Selection
# =============================================================================

if (FDR_MODE) {

  log_step("Part D: FDR calibration across PC levels...")

  files <- sort(list.files(pattern = "^m6A_perm1000.*\\.txt\\.gz$"))
  if (length(files) == 0) stop("No merged permutation files found.")

  numpcs <- data.frame()

  pdf("pval_dists.pdf", useDingbats = FALSE)
  for (f in files) {
    cur <- fread(f, data.table = FALSE)
    q   <- qvalue(cur$V22)
    hist(cur$V22, col = "grey",
         main = paste(f, "\npi1 =", round(1 - q$pi0, 3)))
    pcs <- as.numeric(sub(".*_cov([0-9]+)_merged\\.txt\\.gz$", "\\1", f))
    numpcs <- rbind(numpcs, data.frame(
      no_of_pcs = pcs,
      no_fdr5   = sum(q$qvalues <= 0.05,  na.rm = TRUE),
      no_fdr10  = sum(q$qvalues <= 0.10,  na.rm = TRUE),
      no_fdr25  = sum(q$qvalues <= 0.25,  na.rm = TRUE),
      pi1       = round(1 - q$pi0, 3)
    ))
    message("  PC=", pcs,
            " | FDR5=",  tail(numpcs$no_fdr5,  1),
            " | FDR10=", tail(numpcs$no_fdr10, 1),
            " | pi1=",   tail(numpcs$pi1,      1))
  }
  dev.off()
  message("  p-value distributions → pval_dists.pdf")

  numpcs <- numpcs[order(numpcs$no_of_pcs), ]
  print(numpcs)

  # PC selection plots
  pdf("no_sig_by_pcs_FDRs.pdf", useDingbats = FALSE, height = 7, width = 12)

  y_max <- max(c(numpcs$no_fdr5, numpcs$no_fdr10, numpcs$no_fdr25))
  plot(numpcs$no_of_pcs, numpcs$no_fdr5,
       type = "b", col = "red", pch = 19, ylim = c(0, y_max),
       ylab = "Number of significant hits", xlab = "Number of PCs",
       main = "Significant maQTLs at Different FDR Thresholds")
  lines(numpcs$no_of_pcs, numpcs$no_fdr10, type = "b", col = "blue",  pch = 19)
  lines(numpcs$no_of_pcs, numpcs$no_fdr25, type = "b", col = "green", pch = 19)
  legend("topleft",
         legend = c("FDR 5%", "FDR 10%", "FDR 25%"),
         col    = c("red", "blue", "green"), pch = 19, lty = 1)

  plot(numpcs$no_of_pcs, numpcs$pi1,
       type = "b", pch = 19,
       ylab = expression(pi[1]), xlab = "Number of PCs",
       main = expression("Signal enrichment (" * pi[1] * ") by Number of PCs"))

  dev.off()
  message("  PC selection plots → no_sig_by_pcs_FDRs.pdf")
  message("\n  Review plots, update SELECTED_PC in CONFIG if needed, then re-run --fdr.")

  # Extract significant hits at selected PC level
  selected_file <- paste0("m6A_perm1000_cov", SELECTED_PC, "_merged.txt.gz")
  if (!file.exists(selected_file)) {
    stop("Selected file not found: ", selected_file, "\nCheck SELECTED_PC.")
  }

  if (!file.exists(FDR_SCRIPT)) {
    message("\n  WARNING: ", FDR_SCRIPT, " not found. Run manually:")
    for (fdr in c(0.05, 0.10)) {
      fdr_label  <- ifelse(fdr == 0.05, "FDR5", "FDR10")
      out_prefix <- paste0(OUTPUT_PREFIX, "_", fdr_label, "_", SELECTED_PC, "PCs")
      message("    Rscript ", FDR_SCRIPT, " ", selected_file, " ", fdr, " ", out_prefix)
    }
  } else {
    for (fdr in c(0.05, 0.10)) {
      fdr_label  <- ifelse(fdr == 0.05, "FDR5", "FDR10")
      out_prefix <- paste0(OUTPUT_PREFIX, "_", fdr_label, "_", SELECTED_PC, "PCs")
      message("  Running FDR pass: ", fdr_label)
      system(paste("Rscript", FDR_SCRIPT, selected_file, fdr, out_prefix))
    }
  }

  message("\n  Next: Rscript 01_maQTL_mapping_pipeline.R --nominal")

} # end FDR_MODE

# =============================================================================
# PART E: Submit Nominal Pass
# =============================================================================

if (NOMINAL_MODE) {

  log_step("Part E: Preparing nominal BED and submitting jobs...")

  # Build nominal BED from the permutation BED with two column changes:
  #   pid (col 4): old_pid + "_" + old_gid
  #                e.g. ENST00000360001_1320_GGACT_ENSG00000078808.20
  #   gid (col 5): transcript_id extracted from old_pid (first field before "_")
  #                e.g. ENST00000360001
  # All other columns (coordinates, strand, sample values) stay identical.

  bed_perm <- fread(OUTPUT_BED_GZ, header = TRUE, sep = "\t")
  bed_nom  <- copy(bed_perm)

  old_pid <- bed_nom[[4]]
  old_gid <- bed_nom[[5]]

  # new pid = modification id + "_" + gene id
  bed_nom[[4]] <- paste0(old_pid, "_", old_gid)

  # new gid = transcript id = everything up to and including the first ENST field
  # pid format: ENST00000360001_1320_GGACT  → transcript = ENST00000360001
  bed_nom[[5]] <- sub("_.*", "", old_pid)

  tmp_nom <- sub("\.gz$", "", BED_NOMINAL)
  write.table(bed_nom, tmp_nom, sep = "\t", row.names = FALSE, quote = FALSE)
  system(paste("bgzip -f", tmp_nom))
  system(paste("tabix -p bed", BED_NOMINAL))
  message("  Nominal BED written and indexed: ", BED_NOMINAL)
  message("  Rows: ", nrow(bed_nom), " | pid example: ", bed_nom[[4]][1])
  message("  gid example: ", bed_nom[[5]][1])

  cov_file <- paste0("PC/merged.pcs.", SELECTED_PC)
  system("module load qtltools/1.3")

  for (j in seq_len(N_CHUNKS)) {
    cmd <- paste(
      "QTLtools cis",
      "--vcf", VCF,
      "--bed", BED_NOMINAL,
      "--nominal 1",
      "--cov", cov_file,
      "--normal --grp-best",
      "--chunk", j, N_CHUNKS,
      "--out", paste0("nominals_norm_", j, "_", N_CHUNKS, ".txt")
    )
    submit_lsf(cmd,
               jobname = paste0("nom_", j),
               stdout  = "QTL_nom_norm.out",
               stderr  = "QTL_nom_norm.err")
  }

  message("  Submitted ", N_CHUNKS, " nominal jobs.")
  message("  Wait for completion, then run:")
  message("    Rscript 01_maQTL_mapping_pipeline.R --merge-nominal")

} # end NOMINAL_MODE

# =============================================================================
# PART E (cont.): Merge Nominal Chunks + Gene ID Annotation
# =============================================================================

if (MERGE_NOMINAL) {

  log_step("Part E: Merging nominal chunks and annotating gene IDs...")

  MERGED_FILE    <- "Nominal_RNA_modifications.txt.gz"
  EXTRACTED_FILE <- "extracted_results_all_m6A_modifications_nominal.txt.gz"

  # Merge all chunks
  chunk_files <- sort(list.files(
    pattern = paste0("^nominals_norm_.*_", N_CHUNKS, "\\.txt$")))
  if (length(chunk_files) == 0) stop("No nominal chunk files found.")

  merged <- rbindlist(lapply(chunk_files, fread, header = FALSE))
  fwrite(merged, MERGED_FILE, sep = " ", col.names = FALSE, compress = "gzip")
  message("  Merged ", length(chunk_files), " chunks → ", MERGED_FILE)

  # Extract key columns: mod_id (col 6), SNP (col 10), nominal p-value (col 14)
  system(paste0(
    "zcat ", MERGED_FILE,
    " | cut -d' ' -f6,10,14 | gzip > ", EXTRACTED_FILE
  ))
  message("  Extracted key columns → ", EXTRACTED_FILE)

  # Build gene ID lookup: mod_id → gene_id
  system(paste0("gzip -dc ", OUTPUT_BED_GZ, " | cut -f4,5 > Gene_id_mod_id.txt"))
  message("  Gene ID lookup → Gene_id_mod_id.txt")

  # Join nominals with gene IDs
  system(paste0(
    "zcat ", EXTRACTED_FILE,
    " | grep -F -f <(cut -f1 Gene_id_mod_id.txt)",
    " | sort -k1,1",
    " | join - <(sort -k1,1 Gene_id_mod_id.txt)",
    " > Nominal_gene_ID.txt"
  ))
  n_lines <- as.integer(system("wc -l < Nominal_gene_ID.txt", intern = TRUE))
  message("  Gene-annotated nominals → Nominal_gene_ID.txt (", n_lines, " lines)")

  message("\n  Script 01 complete.")
  message("  Next: Rscript 02_maQTL_replication_QTL.R")

} # end MERGE_NOMINAL
