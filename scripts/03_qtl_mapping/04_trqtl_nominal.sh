#!/usr/bin/env bash
################################################################################
# trQTL Nominal Pass
################################################################################
# Description: Run nominal (all-variant) trQTL analysis to get complete
#              association statistics for all tested transcript-SNP pairs.
#
# Purpose:
#   - Get p-values for ALL transcript-SNP pairs
#   - Required for downstream analyses:
#     * Colocalization with GWAS
#     * Isoform-switching analysis
#     * Trans-QTL discovery
#     * Effect size estimation
#
# Input:  - Transcript expression BED (same as permutation)
#         - Genotype BCF/VCF
#         - Optimal covariate file (determined from permutation)
#
# Output: - Complete nominal association results
#         - Compressed output (larger than eQTL due to more features)
#
# Dependencies: QTLtools, htslib
################################################################################

set -euo pipefail

# Configuration
TRANSCRIPT_BED="/path/to/transcripts_filtered50.bed.gz"
GENOTYPES="/path/to/genotypes.bcf"
OPTIMAL_PCS=0  # UPDATE: Use optimal from permutation step
COV_FILE="/path/to/trQTL/output/pc_optimization/covariates_${OPTIMAL_PCS}PC.txt"
OUTPUT_DIR="/path/to/trQTL/output/nominal"

# QTLtools parameters
CHUNK_SIZE=20
THREADS=1
MEMORY="25000"

# Job submission
PROJECT="acc_yourproject"
QUEUE="premium"

module load qtltools/1.3
module load htslib

mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/chunks"

echo "=========================================="
echo "trQTL Nominal Pass"
echo "=========================================="
echo ""
echo "Transcript Expression: $TRANSCRIPT_BED"
echo "Genotypes: $GENOTYPES"
echo "Covariates (${OPTIMAL_PCS} PCs): $COV_FILE"
echo "Output: $OUTPUT_DIR"
echo ""

################################################################################
# VALIDATION
################################################################################

echo "Validating input files..."

if [[ ! -f "$TRANSCRIPT_BED" ]]; then
    echo "Error: Transcript BED file not found: $TRANSCRIPT_BED"
    exit 1
fi

if [[ ! -f "$GENOTYPES" ]]; then
    echo "Error: Genotype file not found: $GENOTYPES"
    exit 1
fi

if [[ ! -f "$COV_FILE" ]]; then
    echo "Error: Covariate file not found: $COV_FILE"
    echo "Please create it using the permutation script"
    exit 1
fi

echo "All input files found"
echo ""

################################################################################
# RUN NOMINAL PASS IN CHUNKS
################################################################################

echo "Submitting nominal pass jobs..."
echo "----------------------------------------------------------------------"
echo "Note: Nominal pass tests ALL transcript-SNP pairs"
echo "Using --grp-best to group by gene (best transcript per gene)"
echo "This is larger than eQTL due to more transcripts than genes"
echo ""

for chunk in $(seq 1 $CHUNK_SIZE); do

    job_name="trQTL_nominal_${OPTIMAL_PCS}PC_chunk${chunk}"
    out_file="${OUTPUT_DIR}/chunks/trQTL_nominal_${chunk}_${CHUNK_SIZE}.txt"

    echo "QTLtools cis \
        --vcf ${GENOTYPES} \
        --bed ${TRANSCRIPT_BED} \
        --cov ${COV_FILE} \
        --normal \
        --nominal 1 \
        --grp-best \
        --chunk ${chunk} ${CHUNK_SIZE} \
        --out ${out_file}" | \
    bsub -P "$PROJECT" -q "$QUEUE" -n $THREADS -W 20:10 \
        -J "$job_name" \
        -o "$OUTPUT_DIR/logs/${job_name}.out" \
        -e "$OUTPUT_DIR/logs/${job_name}.err" \
        -R "rusage[mem=${MEMORY}]"

    echo "  Submitted chunk ${chunk}/${CHUNK_SIZE}"
done

echo ""
echo "All nominal pass jobs submitted!"
echo "Monitor with: bjobs"
echo ""

read -p "Press Enter when all jobs have completed..."

################################################################################
# MERGE CHUNKS
################################################################################

echo ""
echo "Merging nominal pass results..."
echo "----------------------------------------------------------------------"

# Check if all chunks exist
missing=0
for chunk in $(seq 1 $CHUNK_SIZE); do
    chunk_file="${OUTPUT_DIR}/chunks/trQTL_nominal_${chunk}_${CHUNK_SIZE}.txt"
    if [[ ! -f "$chunk_file" ]]; then
        echo "Warning: Missing chunk file: $chunk_file"
        ((missing++))
    fi
done

if [[ $missing -gt 0 ]]; then
    echo "Error: $missing chunk files are missing"
    echo "Please check job logs in: $OUTPUT_DIR/logs/"
    exit 1
fi

echo "All chunks present. Merging and compressing..."

# Merge all chunks
cat "$OUTPUT_DIR/chunks/trQTL_nominal_"*"_${CHUNK_SIZE}.txt" | \
    gzip -c > "$OUTPUT_DIR/trQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"

# Get file size
file_size=$(du -h "$OUTPUT_DIR/trQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz" | cut -f1)

echo "Merged file created: $OUTPUT_DIR/trQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"
echo "File size: $file_size"
echo ""

# Optional: Remove individual chunks to save space
read -p "Delete individual chunk files to save space? (y/n): " response
if [[ "$response" == "y" ]]; then
    rm "$OUTPUT_DIR/chunks/trQTL_nominal_"*"_${CHUNK_SIZE}.txt"
    echo "Chunk files deleted"
fi

################################################################################
# GENERATE SUMMARY STATISTICS
################################################################################

echo ""
echo "Generating summary statistics..."
echo "----------------------------------------------------------------------"

cat > "$OUTPUT_DIR/summarize_nominal.R" << 'EOF'
#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
nominal_file <- args[1]

cat("Loading nominal results...\n")
nominal <- fread(nominal_file, data.table = FALSE)

cat("\n")
cat("======================================================================\n")
cat("NOMINAL PASS SUMMARY\n")
cat("======================================================================\n")
cat("Total transcript-SNP pairs tested:", nrow(nominal), "\n")
cat("Unique transcripts:", length(unique(nominal$V1)), "\n")
cat("Unique genes:", length(unique(nominal$V2)), "\n")
cat("Unique SNPs:", length(unique(nominal$V8)), "\n")
cat("\n")

# P-value distribution
cat("P-value distribution:\n")
cat("  p < 0.001:", sum(nominal$V12 < 0.001), "\n")
cat("  p < 0.01:", sum(nominal$V12 < 0.01), "\n")
cat("  p < 0.05:", sum(nominal$V12 < 0.05), "\n")
cat("\n")

# Distance to TSS
cat("Distance to TSS:\n")
cat("  Median:", round(median(abs(nominal$V9)), 0), "bp\n")
cat("  Mean:", round(mean(abs(nominal$V9)), 0), "bp\n")
cat("  Max:", max(abs(nominal$V9)), "bp\n")
cat("\n")

# Effect sizes
cat("Effect sizes (beta):\n")
cat("  Median:", round(median(abs(nominal$V13)), 3), "\n")
cat("  Mean:", round(mean(abs(nominal$V13)), 3), "\n")
cat("\n")

cat("======================================================================\n")
EOF

chmod +x "$OUTPUT_DIR/summarize_nominal.R"
Rscript "$OUTPUT_DIR/summarize_nominal.R" "$OUTPUT_DIR/trQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"

echo ""
echo "=========================================="
echo "trQTL Nominal Pass Complete!"
echo "=========================================="
echo ""
echo "Output: $OUTPUT_DIR/trQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"
echo ""
echo "Use this file for:"
echo "  - Colocalization analysis"
echo "  - Isoform-switching detection"
echo "  - Trans-QTL discovery"
echo "  - Comparison with eQTL"
echo ""
