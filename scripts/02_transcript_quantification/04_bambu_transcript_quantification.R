#!/usr/bin/env Rscript
################################################################################
# Bambu Transcript Quantification (Supplementary Notes)
################################################################################
# Description: Alternative transcript quantification using Bambu for comparison
#              with FLAIR (main method). Bambu uses a different isoform
#              discovery and quantification algorithm.
#
# Workflow:
#   1. Discovery: Identify novel isoforms from BAM files
#   2. Quantification: Count reads for all isoforms (known + novel)
#   3. Filtering: Filter protein-coding and lincRNA and apply 50% expression filter
#   4. Classification: Separate annotated vs novel isoforms
#   5. SQANTI QC: Validate discovered isoforms (eventually re-annotate if needed as for FLAIR)
#
# Key Differences from FLAIR:
#   - Uses bayesian framework for isoform assignment
#   - More conservative novel discovery (lower NDR)
#   - Integrated discovery + quantification pipeline
#
# Dependencies: bambu, Rsamtools, GenomeInfoDb, BiocParallel, data.table,
#               ggplot2, gplots, reshape2
################################################################################

suppressPackageStartupMessages({
  library(bambu)
  library(Rsamtools)
  library(GenomeInfoDb)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
  library(gplots)
  library(reshape2)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("BAMBU TRANSCRIPT QUANTIFICATION\n")
cat(rep("=", 70), "\n\n", sep = "")

################################################################################
# SETUP: CONFIGURATION AND ENVIRONMENT
################################################################################

cat("Setting up Bambu environment...\n")
cat(rep("-", 70), "\n", sep = "")

# Fix dbplyr/BiocFileCache compatibility issue
suppressPackageStartupMessages(library(dbplyr))
ns <- asNamespace("dbplyr")
if (exists("db_collect", envir = ns, inherits = FALSE)) {
  unlockBinding("db_collect", ns)
  assign("db_collect", function(con, sql, ...) DBI::dbGetQuery(con, sql), envir = ns)
  lockBinding("db_collect", ns)
}

# Set up output directory and private caches
output_dir <- "/path/to/bambu/output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  ANNOTATIONHUB_CACHE = file.path(output_dir, "AnnotationHub_cache"),
  EXPERIMENTHUB_CACHE = file.path(output_dir, "ExperimentHub_cache")
)
options(ExperimentHub.localHub = TRUE)

# Force single-threaded execution (more stable for large datasets)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1)

cat("Environment configured\n\n")

################################################################################
# STEP 1: DEFINE INPUT FILES
################################################################################

cat("STEP 1: Loading input file paths...\n")
cat(rep("-", 70), "\n", sep = "")

# Reference files
gtf_file <- "/path/to/annotation/gencode.v46lift37.annotation.gtf"
genome_fasta <- "/path/to/reference/hg19.fa"

# BAM files directory
bam_dir <- "/path/to/bam/files"

# All 60 samples
sample_ids <- c(
  "HG00098", "HG00099", "HG00100", "HG00101", "HG00102", "HG00103", "HG00104",
  "HG00110", "HG00111", "HG00112", "HG00115", "HG00116", "HG00117", "HG00118",
  "HG00119", "HG00121", "HG00124", "HG00126", "HG00128", "HG00138", "HG00142",
  "HG00149", "HG00238", "HG00249", "HG00253", "HG00260", "HG00262", "HG00263",
  "HG00265", "NA06985", "NA07037", "NA07051", "NA07357", "NA10847", "NA10851",
  "NA11831", "NA11840", "NA11894", "NA11918", "NA11992", "NA11993", "NA11994",
  "NA12005", "NA12156", "NA12249", "NA12275", "NA12286", "NA12383", "NA12489",
  "NA12750", "NA12760", "NA12761", "NA12762", "NA12763", "NA12776", "NA12812",
  "NA12813", "NA12814", "NA12815", "NA12873"
)

bam_files <- file.path(bam_dir, paste0("sorted_alignment_", sample_ids, ".bam"))

cat("Number of samples:", length(bam_files), "\n")
cat("Reference GTF:", gtf_file, "\n")
cat("Reference genome:", genome_fasta, "\n\n")

################################################################################
# STEP 2: INDEX FILES (IF NEEDED)
################################################################################

cat("STEP 2: Checking and creating indexes...\n")
cat(rep("-", 70), "\n", sep = "")

# Index genome FASTA
if (!file.exists(paste0(genome_fasta, ".fai"))) {
  cat("Indexing genome FASTA...\n")
  Rsamtools::indexFa(genome_fasta)
}

# Index BAM files
invisible(lapply(bam_files, function(bam) {
  if (!file.exists(paste0(bam, ".bai")) && !file.exists(sub("\\.bam$", ".bai", bam))) {
    cat("Indexing", basename(bam), "...\n")
    Rsamtools::indexBam(bam)
  }
}))

cat("All indexes ready\n\n")

################################################################################
# STEP 3: PREPARE ANNOTATIONS
################################################################################

cat("STEP 3: Preparing annotations...\n")
cat(rep("-", 70), "\n", sep = "")

genome <- FaFile(genome_fasta)
annotations <- prepareAnnotations(gtf_file)

# Harmonize chromosome naming (chr1 vs 1)
fa_idx <- Rsamtools::scanFaIndex(genome)
fa_seqs <- as.character(names(fa_idx))

has_chr_fa <- any(startsWith(fa_seqs, "chr"))
has_chr_ann <- any(startsWith(seqlevels(annotations), "chr"))

if (has_chr_fa && !has_chr_ann) {
  cat("Converting annotation to UCSC style (chr prefix)...\n")
  GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"
} else if (!has_chr_fa && has_chr_ann) {
  cat("Converting annotation to NCBI style (no chr prefix)...\n")
  GenomeInfoDb::seqlevelsStyle(annotations) <- "NCBI"
}

# Keep only shared chromosomes
shared_chr <- intersect(seqlevels(annotations), fa_seqs)
annotations <- GenomeInfoDb::keepSeqlevels(annotations, shared_chr, pruning.mode = "coarse")

cat("Annotations prepared:", length(shared_chr), "chromosomes\n\n")

################################################################################
# STEP 4: NOVEL ISOFORM DISCOVERY
################################################################################

cat("STEP 4: Running Bambu isoform discovery...\n")
cat(rep("-", 70), "\n", sep = "")
cat("This step may take 1-3 hours...\n\n")

# Run discovery only (no quantification yet)
# NDR (Novel Discovery Rate) controls stringency
# Higher NDR = more novel isoforms (less stringent)
# Lower NDR = fewer novel isoforms (more stringent)

discovery_result <- bambu(
  reads = bam_files,
  annotations = annotations,
  genome = genome,
  discovery = TRUE,
  quant = FALSE,
  stranded = TRUE,  # RNA-seq is stranded
  ncore = 1,  # Serial for stability
  yieldSize = 1e5,  # Process 100k reads at a time
  lowMemory = TRUE,  # Memory-efficient mode
  rcOutDir = output_dir,
  verbose = TRUE,
  NDR = 0.1  # Novel Discovery Rate (recommended: 0.1)
)

# Save extended annotations (known + novel isoforms)
extended_gtf <- file.path(output_dir, "extended_annotations.gtf")
writeToGTF(discovery_result, file = extended_gtf)

cat("\nDiscovery complete!\n")
cat("Extended GTF saved:", extended_gtf, "\n\n")

################################################################################
# STEP 5: QUANTIFICATION
################################################################################

cat("STEP 5: Quantifying transcript expression...\n")
cat(rep("-", 70), "\n", sep = "")
cat("This step may take 2-4 hours...\n\n")

# Load extended annotations
extended_annotations <- prepareAnnotations(extended_gtf)

# Harmonize chromosome names again
has_chr_ext <- any(startsWith(seqlevels(extended_annotations), "chr"))

if (has_chr_fa && !has_chr_ext) {
  GenomeInfoDb::seqlevelsStyle(extended_annotations) <- "UCSC"
} else if (!has_chr_fa && has_chr_ext) {
  GenomeInfoDb::seqlevelsStyle(extended_annotations) <- "NCBI"
}

# Run quantification
quantification_result <- bambu(
  reads = bam_files,
  annotations = extended_annotations,
  genome = genome,
  discovery = FALSE,  # Use existing annotations
  quant = TRUE,  # Quantify expression
  stranded = TRUE,
  ncore = 1,
  yieldSize = 1e5,
  lowMemory = TRUE,
  rcOutDir = output_dir,
  verbose = TRUE
)

cat("\nQuantification complete!\n\n")

################################################################################
# STEP 6: EXTRACT AND SAVE RESULTS
################################################################################

cat("STEP 6: Extracting expression matrices...\n")
cat(rep("-", 70), "\n", sep = "")

# Extract count and TPM matrices
tx_counts <- assays(quantification_result)$counts
tx_tpm <- assays(quantification_result)$CPM  # Actually TPM in recent Bambu versions

cat("Transcripts quantified:", nrow(tx_counts), "\n")
cat("Samples:", ncol(tx_counts), "\n")

# Get transcript metadata
tx_metadata <- rowRanges(quantification_result)

# Create transcript info table
transcript_info <- data.frame(
  transcript_id = names(tx_metadata),
  gene_id = tx_metadata$GENEID,
  chr = as.character(seqnames(tx_metadata)),
  start = start(tx_metadata),
  end = end(tx_metadata),
  strand = as.character(strand(tx_metadata))
)

# Combine with TPM
tpm_with_info <- cbind(transcript_info, as.data.frame(tx_tpm))

# Save
write.table(tpm_with_info,
            file = file.path(output_dir, "01_bambu_transcript_tpm_all.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("TPM matrix saved\n\n")

################################################################################
# STEP 7: FILTER TO PROTEIN-CODING AND LINCRNA
################################################################################

cat("STEP 7: Filtering to protein-coding and lincRNA genes...\n")
cat(rep("-", 70), "\n", sep = "")

# Load gene types from GTF
gtf_data <- fread(gtf_file, header = FALSE, sep = "\t", skip = 5)
gene_types <- gtf_data[V3 == "gene", .(
  gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
  gene_type = sub('.*gene_type "([^"]+)".*', '\\1', V9)
)]

# Merge with transcript data
tpm_with_type <- merge(
  tpm_with_info,
  gene_types,
  by = "gene_id",
  all.x = TRUE
)

# Filter to protein-coding and lincRNA
tpm_filtered <- tpm_with_type[tpm_with_type$gene_type %in% c("protein_coding", "lincRNA"), ]

cat("Total transcripts:", nrow(tpm_with_info), "\n")
cat("Protein-coding/lincRNA transcripts:", nrow(tpm_filtered), "\n\n")

write.table(tpm_filtered,
            file = file.path(output_dir, "02_transcript_tpm_proteincoding_lincRNA.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# STEP 8: APPLY 50% EXPRESSION FILTER
################################################################################

cat("STEP 8: Applying 50% expression filter...\n")
cat(rep("-", 70), "\n", sep = "")

# Extract annotation columns and TPM columns
setDT(tpm_filtered)
anno_cols <- c("gene_id", "transcript_id", "chr", "start", "end", "strand", "gene_type")
tpm_cols <- setdiff(names(tpm_filtered), anno_cols)

# Calculate proportion of zeros
zero_proportion <- apply(tpm_filtered[, ..tpm_cols], 1, function(x) sum(x == 0) / length(x))

# Keep transcripts with <= 50% zeros
keep <- zero_proportion <= 0.5
tpm_50filter <- tpm_filtered[keep, ]

cat("Before 50% filter:", nrow(tpm_filtered), "transcripts\n")
cat("After 50% filter:", nrow(tpm_50filter), "transcripts\n")
cat("Unique genes:", length(unique(tpm_50filter$gene_id)), "\n\n")

# Reorder columns for BED format
bed_cols <- c("chr", "start", "end", "transcript_id", "gene_id", "strand")
final_cols <- c(bed_cols, tpm_cols)

tpm_50filter_bed <- tpm_50filter[, ..final_cols]
setnames(tpm_50filter_bed, c("#chr", "start", "end", "id", "gid", "strd", tpm_cols))

# Sort by genomic position
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
tpm_50filter_bed$`#chr` <- factor(tpm_50filter_bed$`#chr`, levels = chr_order)
tpm_50filter_sorted <- tpm_50filter_bed[order(`#chr`, start, end)]

write.table(tpm_50filter_sorted,
            file = file.path(output_dir, "03_transcripts_bambu_filtered50.bed"),
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# STEP 9: CLASSIFY ANNOTATED VS NOVEL
################################################################################

cat("STEP 9: Classifying transcripts as annotated vs novel...\n")
cat(rep("-", 70), "\n", sep = "")

# Classify based on transcript ID
# Annotated transcripts have ENST IDs
tpm_50filter_sorted[, is_annotated := grepl("^ENST", id)]
tpm_50filter_sorted[, is_novel := !grepl("^ENST", id)]

# Count per gene
gene_isoform_counts <- tpm_50filter_sorted[, .(
  n_annotated = sum(is_annotated),
  n_novel = sum(is_novel)
), by = gid]

setnames(gene_isoform_counts, "gid", "geneName")

cat("Genes with isoforms:", nrow(gene_isoform_counts), "\n")
cat("  Total annotated isoforms:", sum(gene_isoform_counts$n_annotated), "\n")
cat("  Total novel isoforms:", sum(gene_isoform_counts$n_novel), "\n\n")

write.table(gene_isoform_counts,
            file = file.path(output_dir, "04_gene_isoform_counts.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# STEP 10: GENERATE SUMMARY PLOTS
################################################################################

cat("STEP 10: Generating summary plots...\n")
cat(rep("-", 70), "\n", sep = "")

# Heatmap of annotated vs novel per gene
cat("Creating heatmap...\n")

contingency_table <- table(
  gene_isoform_counts$n_annotated,
  gene_isoform_counts$n_novel
)

color_palette <- colorRampPalette(c("white", "orange", "red"))(n = 299)
color_breaks <- c(
  seq(0, 10, length = 100),
  seq(11, 100, length = 100),
  seq(101, 3000, length = 100)
)

pdf(file.path(output_dir, "05_heatmap_annotated_vs_novel.pdf"), width = 10, height = 8)
heatmap.2(
  contingency_table,
  dendrogram = 'none',
  Rowv = FALSE,
  Colv = FALSE,
  trace = 'none',
  col = color_palette,
  breaks = color_breaks,
  key = FALSE,
  cellnote = ifelse(contingency_table == 0, NA, contingency_table),
  notecol = "black",
  xlab = "Number of Novel Isoforms",
  ylab = "Number of Annotated Isoforms",
  main = "Isoforms per Gene (Bambu, 50% filter)",
  cexRow = 1.2,
  cexCol = 1.2
)
dev.off()

# Histogram: Annotated isoforms per gene
cat("Creating annotated isoforms histogram...\n")

p_annotated <- ggplot(gene_isoform_counts, aes(x = n_annotated)) +
  geom_histogram(breaks = seq(1, max(gene_isoform_counts$n_annotated), by = 1),
                 color = "black", fill = "orange1") +
  ylim(0, max(table(gene_isoform_counts$n_annotated)) * 1.2) +
  labs(
    x = "Number of Annotated Isoforms",
    y = "Number of Genes"
  ) +
  theme_light(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave(file.path(output_dir, "06_annotated_isoforms_per_gene.pdf"),
       p_annotated, height = 8, width = 12)

# Histogram: Novel isoforms per gene
cat("Creating novel isoforms histogram...\n")

p_novel <- ggplot(gene_isoform_counts, aes(x = n_novel)) +
  geom_histogram(breaks = seq(1, max(gene_isoform_counts$n_novel), by = 1),
                 color = "black", fill = "darkgreen") +
  ylim(0, max(table(gene_isoform_counts$n_novel)) * 1.2) +
  labs(
    x = "Number of Novel Isoforms",
    y = "Number of Genes"
  ) +
  theme_light(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave(file.path(output_dir, "07_novel_isoforms_per_gene.pdf"),
       p_novel, height = 8, width = 10)

cat("All plots generated\n\n")

################################################################################
# STEP 11: PREPARE FOR SQANTI VALIDATION (OPTIONAL)
################################################################################

cat("STEP 11: Preparing files for SQANTI3 validation...\n")
cat(rep("-", 70), "\n", sep = "")

# Extract transcript IDs that passed filters
filtered_transcript_ids <- unique(tpm_50filter_sorted$id)

write.table(filtered_transcript_ids,
            file = file.path(output_dir, "transcript_ids_for_sqanti.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Saved", length(filtered_transcript_ids), "transcript IDs\n")
cat("\nTo run SQANTI3:\n")
cat("1. Extract these transcripts from the extended GTF:\n")
cat("   grep -F -f transcript_ids_for_sqanti.txt extended_annotations.gtf > bambu_filtered.gtf\n\n")
cat("2. Run SQANTI3:\n")
cat("   python sqanti3_qc.py bambu_filtered.gtf gencode.gtf hg19.fa\n\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat(rep("=", 70), "\n", sep = "")
cat("BAMBU QUANTIFICATION COMPLETE!\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Summary Statistics:\n")
cat("  Total transcripts discovered:", nrow(tx_counts), "\n")
cat("  Protein-coding/lincRNA:", nrow(tpm_filtered), "\n")
cat("  After 50% filter:", nrow(tpm_50filter_sorted), "\n")
cat("  Genes:", nrow(gene_isoform_counts), "\n")
cat("  Annotated isoforms:", sum(tpm_50filter_sorted$is_annotated), "\n")
cat("  Novel isoforms:", sum(tpm_50filter_sorted$is_novel), "\n\n")

cat("Output Files:\n")
cat("  extended_annotations.gtf - Discovered isoforms\n")
cat("  01_bambu_transcript_tpm_all.txt - All transcripts\n")
cat("  02_transcript_tpm_proteincoding_lincRNA.txt - Filtered by gene type\n")
cat("  03_transcripts_bambu_filtered50.bed - Final for analysis\n")
cat("  04_gene_isoform_counts.txt - Summary per gene\n")
cat("  05-07 PDFs - Summary plots\n\n")

cat("Comparison with FLAIR:\n")
cat("  Use 03_transcripts_bambu_filtered50.bed for comparison\n")
cat("  Compare isoform discovery and quantification\n")
cat(rep("=", 70), "\n", sep = "")
