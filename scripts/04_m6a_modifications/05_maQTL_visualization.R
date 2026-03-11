#!/usr/bin/env Rscript
# =============================================================================
# maQTL Pipeline: Visualization
# =============================================================================
# Description:
#   Produces all main figures:
#
#   1. Residualize modification values
#        Adjust for SEX, GEN, GPC1-3, PC1-3 per modification site via lm().
#
#   2. Violin plots — residualized modification level by genotype
#        One PDF per significant maQTL hit.
#
#   3. Distance plots
#        3a. SNP distance from gene TSS (strand-aware, significant vs all)
#        3b. SNP distance from modification motif (strand-aware)
#
#   4. Transcript count comparison
#        Boxplots: annotated, novel, and total transcripts for
#        mod-genes vs no-mod-genes vs maQTL-genes. Wilcoxon tests.
#
#   5. OAS1 locus LD figure
#        Exon context panel + LD scatter coloured by r² to rs1154970,
#        with eQTL / trQTL / sQTL lead variants highlighted.
#
#   6. OAS1 conditional analysis
#        Effect of rs1154970 on OAS1 m6A before and after conditioning
#        on rs7132797 (eQTL) and rs10774671 (sQTL).
#
# Shell prep required before running (run once):
#   cut -f10 LCLs_m6A_FDR10_3PCs.significant.txt \
#     | sort -u > FDR5_All_SNP.txt
#   bcftools view --include 'ID=@FDR5_All_SNP.txt' genotypes.bcf \
#     -Ov > SNP_all_genotype.txt
#
#   # OAS1 LD file:
#   bcftools view -r chr12:113144767-113544767 -Oz -o oas1_window.vcf.gz genotypes.bcf
#   bcftools index -f oas1_window.vcf.gz
#   plink --vcf oas1_window.vcf.gz --make-bed --out oas1_window
#   plink --bfile oas1_window --r2 --ld-snp rs1154970 \
#     --ld-window 999999 --ld-window-kb 1000 --ld-window-r2 0 \
#     --out oas1_window_LD_rs1154970
#
#   # Dosages for conditional analysis:
#   printf 'rs1154970\nrs7132797\nrs10774671\n' > oas1_3snps.txt
#   bcftools view -i 'ID=@oas1_3snps.txt' genotypes.bcf \
#     | bcftools query -f '%ID[\t%DS]\n' > oas1_3snps.DS.tsv
#   bcftools query -l genotypes.bcf > genotypes.samples.txt
#
# Input:
#   Nanopore_ModificationsRatio_top30_NEW.bed.gz
#   PC/merged.pcs.3
#   LCLs_m6A_FDR10_3PCs.significant.txt
#   m6a_perm1000_cov3_merged.txt.gz
#   SNP_all_genotype.txt
#   Modification_domains_coordinates_ALL.txt
#   Gene_annotated_novel_nonzero_TPM5.txt
#   oas1_window_LD_rs1154970.ld
#   oas1_3snps.DS.tsv / genotypes.samples.txt
#
# Output:
#   Residualized_modifications.txt
#   Dots_SNP_Plots_MaQTL_<GENE>_<SNP>.pdf   (one per hit)
#   MaQTL_all_violin_summary.tsv
#   Distance_TSS_maQTL_All_tpm.pdf
#   Distance_domainmodRNA_maQTL.pdf
#   Annotated_transcripts_boxplot_mod_nomod.pdf
#   Novel_transcripts_boxplot_mod_nomod.pdf
#   All_transcripts_boxplot_mod_nomod.pdf
#   OAS1_window_exons_LD_tagged.pdf
#   OAS1_conditional_analysis_table.txt
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(parallel)
  library(patchwork)
})

# =============================================================================
# CONFIG
# =============================================================================

BED_FILE     <- "Nanopore_ModificationsRatio_top30_NEW.bed.gz"
COV_FILE     <- "PC/merged.pcs.3"
SIG_FILE     <- "LCLs_m6A_FDR10_3PCs.significant.txt"
PERM_FILE    <- "Isoforms_perm1000_cov3_merged.txt.gz"
GENO_FILE    <- "SNP_all_genotype.txt"
MOD_RES_FILE <- "Residualized_modifications.txt"
MOD_COORD    <- "Modification_domains_coordinates_ALL.txt"
GENE_INFO    <- paste0(
  "/sc/arion/projects/bigbrain/Aline_Analysis_Junctions/",
  "FLAIR/5_quantify/Gene_annotated_novel_nonzero_TPM5.txt"
)
GTF_PATH     <- paste0(
  "/sc/arion/projects/bigbrain/Aline_Analysis_Junctions/",
  "FLAIR/Gencode_new/gencode.v46lift37.annotation.gtf"
)
LD_FILE      <- "oas1_window_LD_rs1154970.ld"

N_CORES      <- 6
M6A_POS_OAS1 <- 113346303L   # chr12 genomic position of OAS1 m6A site

META_COLS    <- c("#Chr", "start", "end", "pid", "gid", "strand")

# =============================================================================
# HELPERS
# =============================================================================

log_step <- function(...) message("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

# Parse VCF GT field → allele dosage (0/1/2)
gt_to_dosage <- function(x) {
  x <- gsub(":.*", "", as.character(x))
  x <- gsub("/", "|", x, fixed = TRUE)
  x[x %in% c(".", "./.", ".|.")] <- NA_character_
  x <- gsub("1\\|1",       "2", x)
  x <- gsub("1\\|0|0\\|1", "1", x)
  x <- gsub("0\\|0",       "0", x)
  suppressWarnings(as.numeric(x))
}

# Return sample columns from a VCF data.table row
geno_sample_cols <- function(dt) {
  fixed <- c("#CHROM","CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  setdiff(colnames(dt), fixed)
}

make_violin <- function(df, fill_color = "purple") {
  ggplot(df, aes(x = factor(GT), y = GE)) +
    geom_violin(fill = fill_color, trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_point(size = 3) +
    theme_light(base_size = 18) +
    xlab("Genotype") + ylab("Residualized modification") +
    ggtitle(df$GENE_SNP_STAT[1]) +
    theme(
      axis.text.x  = element_text(size = 28),
      axis.text.y  = element_text(size = 28),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold")
    )
}

make_boxplot_group <- function(data, y_col, ylab_txt, filename) {
  g <- ggplot(data, aes(x = factor(Group), y = .data[[y_col]])) +
    geom_boxplot(aes(fill = factor(Group))) +
    theme_light(base_size = 14) +
    xlab("Gene groups") + ylab(ylab_txt) +
    theme(
      axis.text.x  = element_text(size = 24),
      axis.text.y  = element_text(size = 24),
      axis.title.x = element_text(size = 32, face = "bold"),
      axis.title.y = element_text(size = 32, face = "bold", vjust = 1.8),
      legend.position = "none"
    )
  ggsave(g, filename = filename, height = 8, width = 15, device = "pdf")
  message("  Written: ", filename)
}

# =============================================================================
# STEP 1: Residualize Modifications
# =============================================================================

log_step("Step 1: Residualizing modifications...")

if (file.exists(MOD_RES_FILE)) {
  message("  Found existing file: ", MOD_RES_FILE, " — skipping residualization")
} else {

  Ma6QTL <- as.data.frame(fread(BED_FILE, header = TRUE))
  COV    <- as.data.frame(fread(COV_FILE,  header = TRUE))

  colnames(Ma6QTL)[1] <- "Chr"
  colnames(Ma6QTL)    <- gsub("GM", "NA", colnames(Ma6QTL))
  Ma6QTL[is.na(Ma6QTL)] <- 0

  META_IDX <- seq_along(META_COLS)

  prep_list <- mclapply(seq_len(nrow(Ma6QTL)), function(x) {
    idx <- match(colnames(Ma6QTL)[-META_IDX], colnames(COV))
    data.frame(
      TPM  = as.numeric(Ma6QTL[x, -META_IDX]),
      SEX  = as.factor(COV[1, idx]),
      GEN  = as.factor(COV[2, idx]),
      GPC1 = as.numeric(COV[3, idx]),
      GPC2 = as.numeric(COV[4, idx]),
      GPC3 = as.numeric(COV[5, idx]),
      PC1  = as.numeric(COV[6, idx]),
      PC2  = as.numeric(COV[7, idx]),
      PC3  = as.numeric(COV[8, idx])
    )
  }, mc.cores = N_CORES)
  names(prep_list) <- Ma6QTL$pid

  res_list <- mclapply(prep_list, function(x) {
    residuals(lm(TPM ~ SEX + GEN + GPC1 + GPC2 + GPC3 + PC1 + PC2 + PC3, data = x))
  }, mc.cores = N_CORES)
  names(res_list) <- Ma6QTL$pid

  res_mat <- t(do.call(cbind, res_list))
  out_df  <- cbind(Ma6QTL[, META_IDX], as.data.frame(res_mat, check.names = FALSE))

  write.table(out_df, MOD_RES_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  Written: ", MOD_RES_FILE)
}

# =============================================================================
# STEP 2: Violin Plots — All Significant maQTL Hits
# =============================================================================

log_step("Step 2: Violin plots for all significant maQTL hits...")

if (!file.exists(GENO_FILE)) {
  message("  WARNING: ", GENO_FILE, " not found — skipping violin plots.")
  message("  Generate with:")
  message("    cut -f10 ", SIG_FILE, " | sort -u > FDR5_All_SNP.txt")
  message("    bcftools view --include 'ID=@FDR5_All_SNP.txt' genotypes.bcf -Ov > ", GENO_FILE)
} else {

  sig    <- fread(SIG_FILE,     header = FALSE)
  resmod <- fread(MOD_RES_FILE, header = TRUE)
  geno   <- fread(GENO_FILE,    header = TRUE, check.names = FALSE,
                  skip = function(x) grepl("^##", x))

  summary_rows <- list()

  for (i in seq_len(nrow(sig))) {

    gene_id <- sub("\\..*", "", sig$V1[i])
    snp_id  <- sig$V10[i]

    snp_row <- geno[ID == snp_id]
    if (nrow(snp_row) == 0) { message("  SNP not found: ", snp_id); next }

    g_cols  <- geno_sample_cols(snp_row)
    snp_dos <- sapply(snp_row[1, ..g_cols], gt_to_dosage)

    mod_row <- resmod[grepl(gene_id, pid, fixed = TRUE)][1]
    if (nrow(mod_row) == 0) { message("  Gene not found: ", gene_id); next }

    m_cols   <- setdiff(names(mod_row), META_COLS)
    mod_vals <- setNames(as.numeric(mod_row[1, ..m_cols]), m_cols)

    common <- intersect(names(snp_dos), names(mod_vals))
    if (length(common) < 10) next

    dfp <- data.frame(
      GE = mod_vals[common],
      GT = snp_dos[common]
    )
    dfp <- dfp[is.finite(dfp$GE) & is.finite(dfp$GT), ]

    fit   <- lm(GE ~ GT, data = dfp)
    s     <- summary(fit)
    slope <- unname(fit$coefficients["GT"])
    pval  <- s$coef["GT", "Pr(>|t|)"]

    dfp$GENE_SNP_STAT <- paste0(gene_id, "_", snp_id,
                                "_slope_", round(slope, 2),
                                "_pval_",  signif(pval, 2))

    g <- make_violin(dfp)
    ggsave(g,
           filename = paste0("Dots_SNP_Plots_MaQTL_", gene_id, "_", snp_id, ".pdf"),
           height   = 10, width = 10, device = "pdf", useDingbats = FALSE)

    summary_rows[[length(summary_rows) + 1]] <- data.table(
      Gene_id = gene_id, SNP_ID = snp_id,
      n = nrow(dfp), slope = round(slope, 4), pval = pval
    )
  }

  if (length(summary_rows) > 0) {
    smry <- rbindlist(summary_rows)
    fwrite(smry[order(pval)], "MaQTL_all_violin_summary.tsv", sep = "\t")
    message("  Summary → MaQTL_all_violin_summary.tsv")
  }
}

# =============================================================================
# STEP 3a: Distance from SNP to Gene TSS
# =============================================================================

log_step("Step 3a: SNP distance from gene TSS...")

sig_dt   <- fread(SIG_FILE,  header = FALSE)
perm_dt  <- fread(PERM_FILE, header = FALSE, data.table = TRUE)

# Strand-aware distance: V9 = raw distance, V5 = strand
calc_dtss <- function(dt) {
  d <- dt$V9
  d[dt$V5 == "-"] <- -d[dt$V5 == "-"]
  d
}

tot_sig <- data.frame(
  kbp    = calc_dtss(sig_dt) * 0.001,
  pvalue = sig_dt$V18,
  snp    = sig_dt$V10,
  group  = "significant"
)
tot_all <- data.frame(
  kbp    = calc_dtss(perm_dt) * 0.001,
  pvalue = perm_dt$V18,
  snp    = perm_dt$V10,
  group  = "non-significant"
)
non_sig      <- dplyr::anti_join(tot_all, tot_sig, by = "snp")
combined_tss <- rbind(non_sig, tot_sig)

message("  Significant hits: n=", nrow(tot_sig),
        " | mean |distance| from TSS: ",
        round(mean(abs(tot_sig$kbp), na.rm = TRUE), 2), " kb")

p_tss <- ggplot(combined_tss, aes(x = kbp, y = -log10(pvalue), color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("grey28", "red")) +
  geom_vline(xintercept = 0, color = "red", linewidth = 1, alpha = 0.4) +
  ylab(expression(-log[10](p-value))) +
  xlab("Distance from TSS (kbp)") +
  theme_light(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold")
  )

ggsave(p_tss, filename = "Distance_TSS_maQTL_All_tpm.pdf",
       height = 8, width = 10, device = "pdf")
message("  Written: Distance_TSS_maQTL_All_tpm.pdf")

# =============================================================================
# STEP 3b: Distance from SNP to Modification Motif
# =============================================================================

log_step("Step 3b: SNP distance from modification motif...")

if (!file.exists(MOD_COORD)) {
  message("  WARNING: ", MOD_COORD, " not found — skipping motif distance plot")
} else {

  mod_pos <- fread(MOD_COORD, header = TRUE)

  snp_pos <- data.frame(
    SNP_pos            = sig_dt$V3,
    strand             = sig_dt$V5,
    modification_start = mod_pos$start[match(sig_dt$V6, mod_pos$modification_id)],
    modification_end   = mod_pos$end[  match(sig_dt$V6, mod_pos$modification_id)],
    nominal_pvalue     = sig_dt$V18
  )

  # Strand-aware: positive = SNP is downstream of motif
  snp_pos$distance_bp <- ifelse(
    snp_pos$strand == "+",
    snp_pos$SNP_pos - snp_pos$modification_start,
    snp_pos$modification_end - snp_pos$SNP_pos
  )
  snp_pos$kb <- snp_pos$distance_bp / 1000

  message("  Mean   |distance| to motif: ", round(mean(abs(snp_pos$kb),   na.rm = TRUE), 2), " kb")
  message("  Median |distance| to motif: ", round(median(abs(snp_pos$kb), na.rm = TRUE), 2), " kb")
  message("  Range: ", paste(round(range(abs(snp_pos$kb), na.rm = TRUE), 2), collapse = " – "), " kb")

  p_motif <- ggplot(snp_pos, aes(x = kb, y = -log10(nominal_pvalue))) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, color = "red", linewidth = 1, alpha = 0.4) +
    ylab(expression(-log[10](p-value))) +
    xlab("Distance from modification motif (kb)") +
    theme_light(base_size = 18) +
    theme(
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18),
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold")
    )

  ggsave(p_motif, filename = "Distance_domainmodRNA_maQTL.pdf",
         height = 8, width = 10, device = "pdf")
  message("  Written: Distance_domainmodRNA_maQTL.pdf")

  write.table(snp_pos, "Distance_significant_SNP_modification_motifs.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# =============================================================================
# STEP 4: Transcript Count Comparison — mod genes vs no-mod genes vs maQTL genes
# =============================================================================

log_step("Step 4: Transcript count comparison...")

if (!file.exists(GENE_INFO)) {
  message("  WARNING: GENE_INFO not found — skipping transcript count plots")
} else {

  gene_info <- as.data.frame(fread(GENE_INFO))
  gene_info$geneName <- sub("\\..*", "", gene_info$geneName)

  nano_all_mods <- fread("Nanopore_ModificationsRatio_NAs_filter_17499_clean.bed",
                         header = FALSE)
  mod_genes   <- unique(sub("\\..*", "", nano_all_mods$V5))
  maqtl_genes <- unique(sub("\\..*", "", sig_dt$V1))

  mod   <- merge(data.frame(geneName = mod_genes), gene_info, by = "geneName")
  mod$Group <- "RNA_mod_Genes"
  mod$sum   <- mod$annotated + mod$novel

  no_mod <- gene_info[!gene_info$geneName %in% mod_genes, ]
  no_mod$Group <- "Genes_no_RNA_mod"
  no_mod$sum   <- no_mod$annotated + no_mod$novel

  maqtl_info <- merge(data.frame(geneName = maqtl_genes), gene_info, by = "geneName")
  maqtl_info$Group <- "modQTL_Genes"
  maqtl_info$sum   <- maqtl_info$annotated + maqtl_info$novel

  combined <- rbind(
    mod[,    c("geneName", "annotated", "novel", "Group", "sum")],
    no_mod[, c("geneName", "annotated", "novel", "Group", "sum")]
  )

  make_boxplot_group(combined, "annotated",
                     "Number of annotated transcripts",
                     "Annotated_transcripts_boxplot_mod_nomod.pdf")
  make_boxplot_group(combined, "novel",
                     "Number of novel transcripts",
                     "Novel_transcripts_boxplot_mod_nomod.pdf")
  make_boxplot_group(combined, "sum",
                     "Total number of transcripts",
                     "All_transcripts_boxplot_mod_nomod.pdf")

  message("\n  Wilcoxon tests (mod vs no-mod):")
  for (col in c("annotated", "novel", "sum")) {
    p <- wilcox.test(mod[[col]], no_mod[[col]], alternative = "two.sided")$p.value
    message("  ", col, ": p = ", signif(p, 3))
  }
}

# =============================================================================
# STEP 5: OAS1 Locus LD Figure
# =============================================================================

log_step("Step 5: OAS1 locus LD figure...")

if (!file.exists(LD_FILE)) {
  message("  LD file not found: ", LD_FILE)
  message("  Generate with PLINK (see header notes), then re-run.")
} else {

  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
  })

  ld <- fread(LD_FILE, header = TRUE)
  ld_dt <- ld[, .(pos = BP_B, snp = SNP_B, r2 = R2)]
  ld_dt[, r2_bin := cut(r2,
    breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
    labels = c("0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
    right  = FALSE
  )]
  ld_dt[, x_kb := (pos - M6A_POS_OAS1) / 1000]

  tag_snps <- data.table(
    snp   = c("rs10774671", "rs1154970", "rs7132797"),
    layer = c("sQTL",       "trQTL",     "eQTL")
  )
  ld_dt[, is_tag := snp %in% tag_snps$snp]
  ld_dt <- merge(ld_dt, tag_snps, by = "snp", all.x = TRUE)

  # Exon context from GTF
  gtf_gr <- rtracklayer::import(GTF_PATH)
  ex     <- gtf_gr[gtf_gr$type == "exon" & mcols(gtf_gr)$gene_name == "OAS1"]
  ex     <- ex[seqnames(ex) == "chr12"]
  ex_red <- GenomicRanges::reduce(ex)

  # Find exons flanking the m6A site
  containing <- which(start(ex_red) <= M6A_POS_OAS1 & end(ex_red) >= M6A_POS_OAS1)
  if (length(containing) == 0) {
    i_left  <- max(which(end(ex_red)   < M6A_POS_OAS1))
    i_right <- min(which(start(ex_red) > M6A_POS_OAS1))
    ex2     <- ex_red[c(i_left, i_right)]
  } else {
    i   <- containing[1]
    ex2 <- ex_red[unique(pmax(1, i - 1):pmin(length(ex_red), i + 1))]
  }

  ex_dt <- data.table(start = start(ex2), end = end(ex2))
  ex_dt[, xstart_kb := (start - M6A_POS_OAS1) / 1000]
  ex_dt[, xend_kb   := (end   - M6A_POS_OAS1) / 1000]

  r2_colors <- c("0.0–0.2" = "blue", "0.2–0.4" = "green",
                 "0.4–0.6" = "yellow", "0.6–0.8" = "orange",
                 "0.8–1.0" = "red")

  p_exon <- ggplot() +
    geom_segment(aes(x = min(ld_dt$x_kb), xend = max(ld_dt$x_kb),
                     y = 0, yend = 0),
                 linewidth = 1, color = "grey50") +
    geom_rect(data = ex_dt,
              aes(xmin = xstart_kb, xmax = xend_kb, ymin = -0.25, ymax = 0.25),
              fill = "grey30") +
    geom_vline(xintercept = 0, linetype = 2) +
    annotate("text", x = 0, y = -0.55, label = "m6A site", size = 3) +
    theme_void(base_size = 12) +
    ggtitle("OAS1 exon structure around the m6A site") +
    theme(plot.title = element_text(face = "bold"))

  p_ld <- ggplot(ld_dt, aes(x = x_kb, y = r2, fill = r2_bin)) +
    geom_point(shape = 21, size = 2, color = "black", stroke = 0.15) +
    geom_point(data  = ld_dt[is_tag == TRUE],
               aes(shape = layer), size = 3.5,
               color = "black", fill = "white") +
    geom_text(data  = ld_dt[is_tag == TRUE],
              aes(label = snp), angle = 45, hjust = 0, vjust = 0,
              size = 3, check_overlap = TRUE) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_manual(values = r2_colors, drop = FALSE,
                      name   = expression(r^2 ~ "to anchor")) +
    scale_shape_manual(values = c(eQTL = 16, sQTL = 17, trQTL = 15),
                       na.translate = FALSE) +
    labs(x = "Position relative to m6A site (kb)", y = expression(r^2)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold")) +
    ggtitle("LD to anchor SNP (rs1154970) across OAS1 window")

  fig <- p_exon / p_ld + patchwork::plot_layout(heights = c(0.9, 1.6))
  ggsave("OAS1_window_exons_LD_tagged.pdf", fig,
         width = 7.5, height = 6.5, useDingbats = FALSE)
  message("  Written: OAS1_window_exons_LD_tagged.pdf")
}

# =============================================================================
# STEP 6: OAS1 Conditional Analysis
# =============================================================================

log_step("Step 6: OAS1 conditional analysis...")

cond_inputs <- c("oas1_3snps.DS.tsv", "genotypes.samples.txt")

if (!all(file.exists(cond_inputs))) {
  message("  Dosage files not found — skipping conditional analysis.")
  message("  Generate with (see header notes):")
  message("    printf 'rs1154970\\nrs7132797\\nrs10774671\\n' > oas1_3snps.txt")
  message("    bcftools view -i 'ID=@oas1_3snps.txt' genotypes.bcf \\")
  message("      | bcftools query -f '%ID[\\t%DS]\\n' > oas1_3snps.DS.tsv")
  message("    bcftools query -l genotypes.bcf > genotypes.samples.txt")
} else {

  resmod  <- fread(MOD_RES_FILE, header = TRUE)

  # OAS1 m6A modification values (transcript ENST00000202917)
  oas1_row <- resmod[grepl("ENST00000202917", pid)][1]
  s_cols   <- setdiff(names(oas1_row), META_COLS)
  ph_long  <- data.table(
    Sample = s_cols,
    m6a    = as.numeric(oas1_row[1, ..s_cols])
  )

  # Covariates — reshape wide → long → wide per sample
  cov_wide <- fread(COV_FILE)
  cov_long <- melt(cov_wide, id.vars = names(cov_wide)[1],
                   variable.name = "Sample", value.name = "value")
  setnames(cov_long, names(cov_wide)[1], "Covar")
  cov_mat  <- dcast(cov_long, Sample ~ Covar, value.var = "value")

  df <- merge(ph_long, cov_mat, by = "Sample", all.x = TRUE)

  # Genotype dosages
  geno_raw <- fread("oas1_3snps.DS.tsv", header = FALSE)
  samples  <- fread("genotypes.samples.txt", header = FALSE)[[1]]
  geno_mat <- as.matrix(geno_raw[, -1])
  rownames(geno_mat) <- geno_raw[[1]]
  colnames(geno_mat) <- samples

  geno_df <- data.table(
    Sample     = colnames(geno_mat),
    rs1154970  = as.numeric(geno_mat["rs1154970",  ]),
    rs7132797  = as.numeric(geno_mat["rs7132797",  ]),
    rs10774671 = as.numeric(geno_mat["rs10774671", ])
  )

  df2 <- merge(df, geno_df, by = "Sample")

  covar_candidates <- c("SEX", "GEN", "GPC1", "GPC2", "GPC3",
                        "genes.50percent_1_1_svd_PC1",
                        "genes.50percent_1_1_svd_PC2",
                        "genes.50percent_1_1_svd_PC3")
  covars  <- intersect(covar_candidates, names(df2))
  cov_str <- paste(covars, collapse = " + ")

  get_coef <- function(m, snp) {
    s <- summary(m)$coefficients
    c(beta = s[snp, "Estimate"],
      se   = s[snp, "Std. Error"],
      p    = s[snp, "Pr(>|t|)"])
  }

  models <- list(
    marginal            = lm(as.formula(paste("m6a ~ rs1154970 +", cov_str)), data = df2),
    cond_on_rs7132797   = lm(as.formula(paste("m6a ~ rs1154970 + rs7132797  +", cov_str)), data = df2),
    cond_on_rs10774671  = lm(as.formula(paste("m6a ~ rs1154970 + rs10774671 +", cov_str)), data = df2),
    cond_on_both        = lm(as.formula(paste("m6a ~ rs1154970 + rs7132797 + rs10774671 +", cov_str)), data = df2)
  )

  cond_results <- rbindlist(lapply(names(models), function(nm) {
    r <- get_coef(models[[nm]], "rs1154970")
    data.table(model = nm,
               beta  = round(r["beta"], 4),
               se    = round(r["se"],   4),
               pvalue = r["p"])
  }))

  message("\n  OAS1 rs1154970 conditional analysis:")
  print(cond_results)
  fwrite(cond_results, "OAS1_conditional_analysis_table.txt", sep = "\t")
  message("  Written: OAS1_conditional_analysis_table.txt")
}

message("\nScript 04 complete. All figures generated.")
