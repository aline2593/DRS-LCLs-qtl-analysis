#!/usr/bin/env Rscript
# =============================================================================
# m6A Post-Inference Filtering and Genomic Coordinate Conversion
# =============================================================================
# Description:
#   1. Parse GENCODE v46 annotation to build a transcript-to-gene reference
#   2. Load m6anet merged output (Nanopore_ModificationsRatio.txt)
#   3. Filter modifications:
#        - Remove sites with >50% missing values (NA or probability == 0)
#        - Remove sites with modification probability < 0.9
#        - Keep only protein-coding genes and lncRNAs (Gencode v46)
#   4. Map transcript-relative positions to genomic coordinates
#   5. Output:
#        - Filtered modifications table (.txt)
#        - Sorted genomic coordinates BED file (.bed)
#
# Input files:
#   - filtered_gencode.v46.annotation.txt   : GENCODE annotation (GTF-derived)
#   - Nanopore_ModificationsRatio.txt        : merged m6anet inference output
#
# Output files:
#   - Ref_file_readytouse.txt                  : parsed annotation reference
#   - filtered_m6a_modifications.txt           : filtered modification table
#   - Motifs_position.bed     : genomic coordinates BED file
#
# Usage:
#   Rscript m6a_filter_and_coords.R
# =============================================================================

library(data.table)
library(dplyr)
library(stringr)

# =============================================================================
# CONFIG
# =============================================================================

ANNOTATION_FILE       <- "filtered_gencode.v46.annotation.txt"
MODIFICATIONS_FILE    <- "Nanopore_ModificationsRatio.txt"

REF_OUT               <- "Ref_file_readytouse.txt"
FILTERED_OUT          <- "filtered_m6a_modifications.txt"
BED_OUT               <- "Motifs_position.bed"

PROB_THRESHOLD        <- 0.9    # minimum modification probability to keep
MISSING_THRESHOLD     <- 0.5   # maximum allowed fraction of missing values per site

# Biotypes to retain (as they appear in GENCODE gene_type/gene_biotype field)
KEEP_BIOTYPES <- c("protein_coding", "lncRNA")

# =============================================================================
# STEP 1: Parse GENCODE Annotation
# =============================================================================

message("Step 1: Parsing GENCODE annotation...")

annotation_raw <- fread(ANNOTATION_FILE, header = FALSE, sep = "\t")
message("  Loaded annotation: ", nrow(annotation_raw), " rows x ", ncol(annotation_raw), " cols")

# Split the semicolon-delimited attribute field (V9)
split_cols  <- strsplit(as.character(annotation_raw[["V9"]]), ";")
max_fields  <- max(sapply(split_cols, length))
col_names   <- paste0("V9_", seq_len(max_fields))

split_df <- as.data.frame(
  do.call(rbind, lapply(split_cols, function(x) c(x, rep(NA, max_fields - length(x)))))
)
colnames(split_df) <- col_names

annotation <- cbind(annotation_raw[, -"V9"], split_df)

# Keep columns 1:9 (chr, source, feature, start, end, score, strand, frame,
# gene_id, transcript_id)
annotation <- annotation[, 1:9]
colnames(annotation)[8:9] <- c("gene_id_raw", "transcript_id_raw")

# Clean gene_id and transcript_id: strip key name, version suffix, whitespace, quotes
clean_id <- function(x, key) {
  x <- gsub(key, "", x)          # remove key name
  x <- sub("\\..*", "", x)       # remove version suffix (.1, .2, ...)
  x <- trimws(x)                  # strip whitespace
  x <- sub('^"', "", x)           # strip leading quote if present
  x <- sub('"$', "", x)           # strip trailing quote if present
  x
}

annotation$gene_id       <- clean_id(annotation$gene_id_raw,       "gene_id ")
annotation$transcript_id <- clean_id(annotation$transcript_id_raw, "transcript_id ")

ref <- annotation[, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "gene_id", "transcript_id")]
colnames(ref)[1:7] <- c("chr", "source", "feature", "start", "end", "score", "strand")

write.table(ref, REF_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
message("  Reference file written: ", REF_OUT)

# =============================================================================
# STEP 2: Load m6anet Modifications
# =============================================================================

message("Step 2: Loading m6anet modifications...")

mod_raw <- fread(MODIFICATIONS_FILE, header = TRUE, sep = "\t")
message("  Loaded modifications: ", nrow(mod_raw), " rows x ", ncol(mod_raw), " cols")

# Keep core columns (first 6: transcript_id_pos, chr, start, id, probability cols)
mod <- mod_raw[, 1:6]

# Rename column 4 to 'id' (transcript_id_position field e.g. ENST00000_123)
colnames(mod)[4] <- "id"

# Split id into transcript prefix and position suffix
split_id        <- str_split_fixed(mod$id, "_", n = 2)
mod$transcript  <- split_id[, 1]   # e.g. ENST00000123456
mod$position    <- as.numeric(split_id[, 2])  # position within transcript

message("  Parsed transcript IDs and positions.")

# =============================================================================
# STEP 3: Filter Modifications
# =============================================================================

message("Step 3: Filtering modifications...")
n_start <- nrow(mod)

# --- 3a. Identify probability columns (all numeric cols except position/coords)
# Assumes probability columns are columns 5 and 6 onward in the original file,
# or any numeric column that isn't a coordinate. Adjust if your file differs.
prob_cols <- colnames(mod_raw)[5:ncol(mod_raw)]
prob_cols <- prob_cols[prob_cols %in% colnames(mod)]

# Rebuild with all probability columns for filtering
mod_full <- mod_raw[, 1:6]
colnames(mod_full)[4] <- "id"
split_id2 <- str_split_fixed(mod_full$id, "_", n = 2)
mod_full$transcript <- split_id2[, 1]
mod_full$position   <- as.numeric(split_id2[, 2])

# --- 3b. Remove sites where >50% of samples have NA or probability == 0
# (both count as missing: NA = site not detected; 0 = no modification signal)
prob_data <- mod_raw[, 5:ncol(mod_raw), drop = FALSE]
n_samples <- ncol(prob_data)

missing_per_site <- rowSums(is.na(prob_data) | prob_data == 0) / n_samples
pass_missing <- missing_per_site <= MISSING_THRESHOLD

mod_full <- mod_full[pass_missing, ]
message("  After missing value filter (>", MISSING_THRESHOLD * 100, "% missing): ",
        nrow(mod_full), " sites (removed ", n_start - nrow(mod_full), ")")

# --- 3c. Remove sites with modification probability < threshold
# Column 5 in the original file is the site-level probability (n_mod_reads / n_reads)
# and column 6 is the probability_modified from m6anet. Adjust col name if needed.
prob_col_name <- colnames(mod_raw)[6]   # typically "probability_modified"
mod_full <- mod_full[mod_full[[prob_col_name]] >= PROB_THRESHOLD, ]
message("  After probability filter (>=", PROB_THRESHOLD, "): ",
        nrow(mod_full), " sites (removed ", n_start - nrow(mod_full), " total so far)")

# =============================================================================
# STEP 4: Annotate with Gene Info and Filter by Biotype
# =============================================================================

message("Step 4: Annotating modifications with gene/biotype info...")

# Merge modifications with reference on transcript_id
merged <- merge(mod_full, ref, by.x = "transcript", by.y = "transcript_id", all.x = FALSE)
message("  After merge with annotation: ", nrow(merged), " rows")

# Filter to protein-coding and lncRNA biotypes
# The 'source' field in a filtered GTF often carries the biotype for transcript rows
# If your annotation has a dedicated gene_biotype column, use that instead
merged_filtered <- merged[merged$source %in% KEEP_BIOTYPES, ]
message("  After biotype filter (", paste(KEEP_BIOTYPES, collapse = ", "), "): ",
        nrow(merged_filtered), " modifications on ",
        length(unique(merged_filtered$gene_id)), " genes")

# Write filtered table
out_cols <- c("transcript", "position", colnames(mod_raw)[5:6], "gene_id", "chr", "start", "end", "strand", "source")
out_cols <- out_cols[out_cols %in% colnames(merged_filtered)]
write.table(merged_filtered[, out_cols, with = FALSE],
            FILTERED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
message("  Filtered modifications written: ", FILTERED_OUT)

# =============================================================================
# STEP 5: Map Transcript Positions to Genomic Coordinates
# =============================================================================

message("Step 5: Computing genomic coordinates...")

# modification_position = transcript annotation start + position within transcript
merged_filtered$modification_start <- as.numeric(merged_filtered$start) +
                                       as.numeric(merged_filtered$position)
merged_filtered$modification_end   <- merged_filtered$modification_start + 5

# Build BED format: chr, start, end, id, probability, strand
bed <- merged_filtered[, c("chr", "modification_start", "modification_end",
                            "transcript", "gene_id", "strand")]
colnames(bed) <- c("#chr", "start", "end", "transcript_id", "gene_id", "strand")

# Sort by chromosome and position (standard genomic order)
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
bed_sorted <- bed %>%
  mutate(`#chr` = factor(`#chr`, levels = chr_order)) %>%
  arrange(`#chr`, start, end) %>%
  mutate(`#chr` = as.character(`#chr`))

write.table(bed_sorted, BED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
message("  BED file written: ", BED_OUT)

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Summary ===")
message("  Input modifications:          ", n_start)
message("  After missing value filter:   ", sum(pass_missing))
message("  After probability filter:     ", nrow(mod_full))
message("  After biotype filter:         ", nrow(merged_filtered))
message("  Unique genes retained:        ", length(unique(merged_filtered$gene_id)))
message("  Output files:")
message("    - ", REF_OUT)
message("    - ", FILTERED_OUT)
message("    - ", BED_OUT)
message("\nDone.")
