#!/usr/bin/env bash
################################################################################
# StringTie Transcript Quantification (Supplementary Notes)
################################################################################
# Description: Alternative transcript quantification using StringTie with
#              transcript trimming option for comparison with FLAIR and Bambu.
#
# Workflow:
#   1. Per-sample assembly: Discover isoforms in each sample
#   2. Merge: Combine all sample GTFs into unified annotation
#   3. Quantification: Re-quantify all samples using merged annotation
#   4. Filtering: Extract counts/TPM and apply filters
#   5. SQANTI validation: Validate discovered isoforms
#
# Key Parameters (Trimmed Version - Final):
#   -c 1.5    : Minimum coverage threshold (TPM)
#   -t        : Trim predicted transcripts (remove artifacts at ends)
#   -g 50     : Minimum gap between transcripts (bp)
#   -L        : Label prefix for novel transcripts (MSTRG)
#
# Dependencies: StringTie, samtools, R (rtracklayer, data.table, ggplot2)
################################################################################

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

################################################################################
# CONFIGURATION
################################################################################

# Input/output paths
BAM_DIR="/path/to/bam/files"
GTF_FILE="/path/to/annotation/gencode.v46lift37.annotation.gtf"
GENOME_FA="/path/to/reference/hg19.fa"
OUTPUT_DIR="/path/to/stringtie/output"
POLYA_MOTIFS="/path/to/mouse_and_human.polyA_motif.txt"

# StringTie parameters (TRIMMED VERSION - FINAL)
THREADS=8
MIN_COVERAGE=1.5  # Minimum TPM for transcript
MIN_GAP=50        # Minimum gap between transcripts (bp)
USE_TRIMMING=true # Trim transcript ends (removes artifacts)

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

echo "Configuration:"
echo "  BAM directory: $BAM_DIR"
echo "  Reference GTF: $GTF_FILE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Threads: $THREADS"
echo "  Min coverage: $MIN_COVERAGE"
echo "  Min gap: $MIN_GAP"
echo "  Trimming: $USE_TRIMMING"
echo ""

################################################################################
# STEP 1: PER-SAMPLE TRANSCRIPT ASSEMBLY
################################################################################

sample_count=0
for bam in "$BAM_DIR"/*.bam; do
    sample_name=$(basename "$bam" .bam)

    # Skip if output already exists
    if [ -f "$OUTPUT_DIR/${sample_name}.gtf" ]; then
        echo "  Skipping $sample_name (already processed)"
        continue
    fi

    echo "  Processing: $sample_name"

    # Run StringTie assembly with trimming
    stringtie "$bam" \
        -G "$GTF_FILE" \
        -o "$OUTPUT_DIR/${sample_name}.gtf" \
        -p "$THREADS" \
        -c "$MIN_COVERAGE" \
        -t \
        -L \
        -g "$MIN_GAP" \
        -A "$OUTPUT_DIR/${sample_name}_gene_abundance.tsv" \
        2>> "$OUTPUT_DIR/logs/${sample_name}_assembly.log"

    ((sample_count++))
done

echo ""
echo "Assembly complete for $sample_count samples"
echo ""

################################################################################
# STEP 2: MERGE ALL SAMPLE GTFs
################################################################################

mkdir -p "$OUTPUT_DIR/merged"

# Create list of sample GTFs
find "$OUTPUT_DIR" -maxdepth 1 -name "*.gtf" \
    ! -name "merged_transcripts.gtf" \
    > "$OUTPUT_DIR/gtf_list.txt"

n_gtfs=$(wc -l < "$OUTPUT_DIR/gtf_list.txt")
echo "  Found $n_gtfs sample GTFs to merge"

# Merge all GTFs
stringtie --merge \
    -G "$GTF_FILE" \
    -o "$OUTPUT_DIR/merged/merged_transcripts.gtf" \
    "$OUTPUT_DIR/gtf_list.txt" \
    2>> "$OUTPUT_DIR/logs/merge.log"

echo "  Merged GTF saved: $OUTPUT_DIR/merged/merged_transcripts.gtf"
echo ""

################################################################################
# STEP 3: RE-QUANTIFY ALL SAMPLES WITH MERGED ANNOTATION
################################################################################

mkdir -p "$OUTPUT_DIR/ballgown"
mkdir -p "$OUTPUT_DIR/logs"

quant_count=0

for bam in "$BAM_DIR"/*.bam; do
    sample_name=$(basename "$bam" .bam)

    # Skip if already quantified
    if [ -f "$OUTPUT_DIR/ballgown/${sample_name}/${sample_name}_quant.gtf" ]; then
        echo "  Skipping $sample_name (already quantified)"
        continue
    fi

    echo "  Quantifying: $sample_name"

    # Create per-sample ballgown directory
    # Required for -B flag — each sample needs its own directory
    # otherwise ballgown tables overwrite each other
    mkdir -p "$OUTPUT_DIR/ballgown/${sample_name}"

    stringtie "$bam" \
        -G "$OUTPUT_DIR/merged/merged_transcripts.gtf" \
        -o "$OUTPUT_DIR/ballgown/${sample_name}/${sample_name}_quant.gtf" \
        -e \
        -B \
        -L \
        -p "$THREADS" \
        2>> "$OUTPUT_DIR/logs/${sample_name}_quant.log"

    quant_count=$((quant_count + 1))   # safe with set -e

    echo "  Done: $sample_name"
done

echo ""
echo "Quantification complete for $quant_count samples"
echo ""

################################################################################
# STEP 4: EXTRACT TRANSCRIPT COUNTS AND TPM
################################################################################

echo "Creating sample list for prepDE..."

# Build sample list pointing to per-sample ballgown subdirectories
ls "$OUTPUT_DIR"/ballgown/*/*.gtf | while read file; do
    sample=$(basename "$(dirname "$file")")
    echo -e "${sample}\t${file}"
done > "$OUTPUT_DIR/sample_list.tsv"

echo "Sample list created: $OUTPUT_DIR/sample_list.tsv"
echo "Samples found: $(wc -l < $OUTPUT_DIR/sample_list.tsv)"

# Check if prepDE.py3 exists
PREPDE_SCRIPT="/sc/arion/projects/bigbrain/Aline_Analysis_Junctions/FLAIR/Stringtie/Trimm/prepDE.py3"

if [ ! -f "$PREPDE_SCRIPT" ]; then
    echo "Error: prepDE.py3 not found at $PREPDE_SCRIPT"
    echo "Please download from: https://github.com/gpertea/stringtie/blob/master/prepDE.py3"
    exit 1
fi

# Compute median read length from first BAM file
# This is used by prepDE to convert coverage to read counts
# For long reads this should reflect actual median read length 
FIRST_BAM=$(ls "$BAM_DIR"/*.bam | head -1)
echo "Computing median read length from: $(basename $FIRST_BAM)"
MEDIAN_READ_LENGTH=$(samtools view "$FIRST_BAM" \
    | awk '{print length($10)}' \
    | sort -n \
    | awk 'BEGIN{c=0} {a[c++]=$1} END{print a[int(c/2)]}')
echo "Median read length: $MEDIAN_READ_LENGTH bp"

# Run prepDE
echo "Extracting count matrices..."
python "$PREPDE_SCRIPT" \
    -i "$OUTPUT_DIR/sample_list.tsv" \
    -g "$OUTPUT_DIR/gene_count_matrix.csv" \
    -t "$OUTPUT_DIR/transcript_count_matrix.csv" \
    -l "$MEDIAN_READ_LENGTH" \
    2>> "$OUTPUT_DIR/logs/prepDE.log"

echo "Count matrices created:"
echo "  Gene:       $OUTPUT_DIR/gene_count_matrix.csv"
echo "  Transcript: $OUTPUT_DIR/transcript_count_matrix.csv"
echo ""

################################################################################
# STEP 5: PROCESS IN R - CONVERT TO TPM AND FILTER
################################################################################

# Create R script
cat > "$OUTPUT_DIR/process_stringtie.R" << 'EOF'
#!/usr/bin/env Rscript
################################################################################
# StringTie Post-Processing
################################################################################

suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gplots)
  library(reshape2)
})

setwd(Sys.getenv("WORK_DIR"))

cat("\nProcessing StringTie results...\n")
cat("======================================================================\n\n")

################################################################################
# Load merged GTF and extract transcript info
################################################################################

cat("Loading merged GTF...\n")

gtf_file <- "merged/merged_transcripts.gtf"
gtf_data <- import(gtf_file, format = "gtf")

# Convert to data frame and filter for transcripts
gtf_df <- as.data.frame(gtf_data)
transcript_df <- gtf_df %>%
  filter(type == "transcript") %>%
  select(seqnames, start, end, strand, gene_id, transcript_id, gene_name)

cat("Total transcripts in merged GTF:", nrow(transcript_df), "\n")

# Classify as annotated or novel
transcript_df <- transcript_df %>%
  mutate(type = ifelse(grepl("^ENST", transcript_id), "Annotated", "Novel"))

cat("  Annotated:", sum(transcript_df$type == "Annotated"), "\n")
cat("  Novel:", sum(transcript_df$type == "Novel"), "\n\n")

################################################################################
# Load count matrix and convert to TPM
################################################################################

cat("Loading count matrix...\n")

counts <- fread("transcript_count_matrix.csv", data.table = FALSE)
colnames(counts)[1] <- "transcript_id"

cat("Transcripts with counts:", nrow(counts), "\n")
cat("Samples:", ncol(counts) - 1, "\n\n")

# Merge with transcript info
counts_with_info <- merge(transcript_df, counts, by = "transcript_id")

# Calculate transcript lengths
counts_with_info$length <- abs(counts_with_info$end - counts_with_info$start)

# Convert counts to TPM
cat("Converting counts to TPM...\n")

sample_cols <- grep("^sorted_alignment_", names(counts_with_info), value = TRUE)
count_matrix <- as.matrix(counts_with_info[, sample_cols])

# TPM = (counts / length) / sum(counts / length) * 1e6
rate <- count_matrix / (counts_with_info$length / 1000)  # RPK
tpm <- t(t(rate) / colSums(rate)) * 1e6

# Combine with metadata
tpm_df <- cbind(counts_with_info[, c("seqnames", "start", "end", "strand",
                                     "gene_id", "transcript_id", "gene_name", "type")],
                as.data.frame(tpm))

cat("TPM conversion complete\n\n")

# Save all transcripts with TPM
fwrite(tpm_df, "01_all_transcripts_tpm.txt", sep = "\t", quote = FALSE)

################################################################################
# Filter to protein-coding and lincRNA genes
################################################################################

cat("Filtering to protein-coding and lincRNA genes...\n")

# Load gene types from GTF
gtf_genes <- gtf_df %>%
  filter(type == "gene") %>%
  select(gene_id, gene_type) %>%
  distinct()

# Merge and filter
tpm_with_type <- merge(tpm_df, gtf_genes, by = "gene_id", all.x = TRUE)
tpm_pc_lnc <- tpm_with_type %>%
  filter(gene_type %in% c("protein_coding", "lincRNA"))

cat("Protein-coding/lincRNA transcripts:", nrow(tpm_pc_lnc), "\n\n")

fwrite(tpm_pc_lnc, "02_transcripts_proteincoding_lincRNA.txt", sep = "\t", quote = FALSE)

################################################################################
# Apply 50% expression filter
################################################################################

cat("Applying 50% expression filter...\n")

setDT(tpm_pc_lnc)

# Calculate proportion of zeros
sample_cols_final <- grep("^sorted_alignment_", names(tpm_pc_lnc), value = TRUE)
tpm_matrix <- as.matrix(tpm_pc_lnc[, ..sample_cols_final])

zero_proportion <- apply(tpm_matrix, 1, function(x) sum(x == 0) / length(x))

# Keep transcripts expressed in >= 50% of samples
keep <- zero_proportion <= 0.5
tpm_filtered <- tpm_pc_lnc[keep, ]

cat("Before 50% filter:", nrow(tpm_pc_lnc), "\n")
cat("After 50% filter:", nrow(tpm_filtered), "\n")
cat("Unique genes:", length(unique(tpm_filtered$gene_id)), "\n\n")

# Format as BED
bed_cols <- c("seqnames", "start", "end", "transcript_id", "gene_id", "strand")
final_cols <- c(bed_cols, sample_cols_final)

bed_df <- tpm_filtered[, ..final_cols]
setnames(bed_df, c("#chr", "start", "end", "id", "gid", "strd", sample_cols_final))

# Sort by genomic position
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
bed_df$`#chr` <- factor(bed_df$`#chr`, levels = chr_order)
bed_sorted <- bed_df[order(`#chr`, start, end)]

fwrite(bed_sorted, "03_transcripts_filtered50.bed", sep = "\t", quote = FALSE)

cat("BED file saved for QTL analysis\n\n")

################################################################################
# Count isoforms per gene
################################################################################

cat("Counting isoforms per gene...\n")

# Classify by type
bed_sorted[, is_annotated := grepl("^ENST", id)]
bed_sorted[, is_novel := grepl("^MSTRG", id)]

gene_counts <- bed_sorted[, .(
  n_annotated = sum(is_annotated),
  n_novel = sum(is_novel)
), by = gid]

setnames(gene_counts, "gid", "geneName")

cat("Genes with isoforms:", nrow(gene_counts), "\n")
cat("  Total annotated:", sum(gene_counts$n_annotated), "\n")
cat("  Total novel:", sum(gene_counts$n_novel), "\n\n")

fwrite(gene_counts, "04_gene_isoform_counts.txt", sep = "\t", quote = FALSE)

################################################################################
# Generate summary plots
################################################################################

cat("Generating summary plots...\n")

# Heatmap
contingency <- table(gene_counts$n_annotated, gene_counts$n_novel)

color_palette <- colorRampPalette(c("white", "orange", "red"))(n = 299)
color_breaks <- c(
  seq(0, 10, length = 100),
  seq(11, 100, length = 100),
  seq(101, 3000, length = 100)
)

pdf("05_heatmap_annotated_vs_novel.pdf", width = 10, height = 8)
heatmap.2(
  contingency,
  dendrogram = 'none',
  Rowv = FALSE,
  Colv = FALSE,
  trace = 'none',
  col = color_palette,
  breaks = color_breaks,
  key = FALSE,
  cellnote = ifelse(contingency == 0, NA, contingency),
  notecol = "black",
  xlab = "Number of Novel Isoforms",
  ylab = "Number of Annotated Isoforms",
  main = "StringTie Isoforms per Gene (50% filter)",
  cexRow = 1.2,
  cexCol = 1.2
)
dev.off()

# Histogram: Annotated
p_anno <- ggplot(gene_counts, aes(x = n_annotated)) +
  geom_histogram(breaks = seq(1, max(gene_counts$n_annotated), by = 1),
                 color = "black", fill = "orange1") +
  labs(x = "Number of Annotated Isoforms", y = "Number of Genes") +
  theme_light(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave("06_annotated_isoforms_per_gene.pdf", p_anno, height = 8, width = 12)

# Histogram: Novel
p_novel <- ggplot(gene_counts, aes(x = n_novel)) +
  geom_histogram(breaks = seq(1, max(gene_counts$n_novel), by = 1),
                 color = "black", fill = "darkgreen") +
  labs(x = "Number of Novel Isoforms", y = "Number of Genes") +
  theme_light(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

ggsave("07_novel_isoforms_per_gene.pdf", p_novel, height = 8, width = 10)

cat("Plots generated\n\n")

################################################################################
# Prepare for SQANTI
################################################################################

cat("Preparing transcript list for SQANTI...\n")

writeLines(unique(bed_sorted$id), "transcript_ids_for_sqanti.txt")

cat("Saved", length(unique(bed_sorted$id)), "transcript IDs\n\n")

cat("======================================================================\n")
cat("Processing complete!\n")
cat("======================================================================\n")
EOF

# Run R script
echo "  Running R processing script..."
WORK_DIR="$OUTPUT_DIR" Rscript "$OUTPUT_DIR/process_stringtie.R"

echo ""

################################################################################
# STEP 6: PREPARE FOR SQANTI VALIDATION
################################################################################

echo "  Extracting filtered transcripts from GTF..."

# Filter GTF to keep only transcripts that passed filters
grep -F -f "$OUTPUT_DIR/transcript_ids_for_sqanti.txt" \
    "$OUTPUT_DIR/merged/merged_transcripts.gtf" \
    > "$OUTPUT_DIR/filtered_transcripts_for_sqanti.gtf"

echo "  Filtered GTF created: $OUTPUT_DIR/filtered_transcripts_for_sqanti.gtf"
echo ""

################################################################################
# STEP 7: RUN SQANTI (OPTIONAL)
################################################################################

if command -v python &> /dev/null && [ -d "/path/to/SQANTI3" ]; then

    module load sqanti3 || true

    python /path/to/SQANTI3/sqanti3_qc.py \
        "$OUTPUT_DIR/filtered_transcripts_for_sqanti.gtf" \
        "$GTF_FILE" \
        "$GENOME_FA" \
        --force_id_ignore \
        --polyA_motif_list "$POLYA_MOTIFS" \
        -d "$OUTPUT_DIR/sqanti_output" \
        --report both

    echo "  SQANTI results: $OUTPUT_DIR/sqanti_output"
    echo ""
else
    echo "Skipping SQANTI (not configured)"
    echo "To run manually:"
    echo "  python sqanti3_qc.py filtered_transcripts_for_sqanti.gtf \\"
    echo "    gencode.gtf hg19.fa --force_id_ignore --report both"
    echo ""
fi

################################################################################
# FINAL SUMMARY
################################################################################

echo ""
echo "Output Files:"
echo "  merged/merged_transcripts.gtf - Merged annotation"
echo "  01_all_transcripts_tpm.txt - All transcripts with TPM"
echo "  02_transcripts_proteincoding_lincRNA.txt - Filtered by gene type"
echo "  03_transcripts_filtered50.bed - Final for QTL analysis"
echo "  04_gene_isoform_counts.txt - Isoforms per gene"
echo "  05-07 PDFs - Summary plots"
echo "  sqanti_output/ - SQANTI validation results"
echo ""
echo "Next Steps:"
echo "  1. Compress BED file:"
echo "     bgzip 03_transcripts_filtered50.bed"
echo "     tabix -p bed 03_transcripts_filtered50.bed.gz"
echo ""
echo "  2. Compare with FLAIR and Bambu results"
echo ""
