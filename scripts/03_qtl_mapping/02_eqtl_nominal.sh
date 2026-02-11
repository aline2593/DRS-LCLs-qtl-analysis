#!/usr/bin/env bash
################################################################################
# eQTL Nominal Pass
################################################################################
# Description: Run nominal (all-variant) eQTL analysis to get complete
#              association statistics for all tested gene-SNP pairs.
#
# Purpose:
#   - Get p-values for ALL gene-SNP pairs (not just best per gene)
#   - Required for downstream analyses:
#     * Fine-mapping
#     * Colocalization with GWAS
#     * Mendelian randomization
#     * Effect size estimation
#
# Input:  - Gene expression BED (same as permutation)
#         - Genotype BCF/VCF
#         - Optimal covariate file (determined from permutation)
#
# Output: - Complete nominal association results for all gene-SNP pairs
#         - Compressed output (can be very large)
#
# Dependencies: QTLtools, htslib
################################################################################

set -euo pipefail

# Configuration
EXPRESSION_BED="/path/to/gene_expression_filtered50.bed.gz"
GENOTYPES="/path/to/genotypes.bcf"
OPTIMAL_PCS=3  # UPDATE: Use optimal from permutation step
COV_FILE="/path/to/eqtl/output/pc_optimization/covariates_${OPTIMAL_PCS}PC.txt"
OUTPUT_DIR="/path/to/eqtl/output/nominal"

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
echo "eQTL Nominal Pass"
echo "=========================================="
echo ""
echo "Expression: $EXPRESSION_BED"
echo "Genotypes: $GENOTYPES"
echo "Covariates (${OPTIMAL_PCS} PCs): $COV_FILE"
echo "Output: $OUTPUT_DIR"
echo ""

################################################################################
# VALIDATION
################################################################################

echo "Validating input files..."

if [[ ! -f "$EXPRESSION_BED" ]]; then
    echo "Error: Expression file not found: $EXPRESSION_BED"
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
echo "Note: Nominal pass tests ALL gene-SNP pairs"
echo "This generates large output files and takes longer than permutation"
echo ""

for chunk in $(seq 1 $CHUNK_SIZE); do

    job_name="eQTL_nominal_${OPTIMAL_PCS}PC_chunk${chunk}"
    out_file="${OUTPUT_DIR}/chunks/eQTL_nominal_${chunk}_${CHUNK_SIZE}.txt"

    echo "QTLtools cis \
        --vcf ${GENOTYPES} \
        --bed ${EXPRESSION_BED} \
        --cov ${COV_FILE} \
        --normal \
        --nominal 1 \
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
    chunk_file="${OUTPUT_DIR}/chunks/eQTL_nominal_${chunk}_${CHUNK_SIZE}.txt"
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
cat "$OUTPUT_DIR/chunks/eQTL_nominal_"*"_${CHUNK_SIZE}.txt" | \
    gzip -c > "$OUTPUT_DIR/eQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"

# Get file size
file_size=$(du -h "$OUTPUT_DIR/eQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz" | cut -f1)

echo "Merged file created: $OUTPUT_DIR/eQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"
echo "File size: $file_size"
echo ""

# Optional: Remove individual chunks to save space
read -p "Delete individual chunk files to save space? (y/n): " response
if [[ "$response" == "y" ]]; then
    rm "$OUTPUT_DIR/chunks/eQTL_nominal_"*"_${CHUNK_SIZE}.txt"
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
cat("Total gene-SNP pairs tested:", nrow(nominal), "\n")
cat("Unique genes:", length(unique(nominal$V1)), "\n")
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

cat("======================================================================\n")
EOF

chmod +x "$OUTPUT_DIR/summarize_nominal.R"
Rscript "$OUTPUT_DIR/summarize_nominal.R" "$OUTPUT_DIR/eQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"

echo ""
echo "=========================================="
echo "eQTL Nominal Pass Complete!"
echo "=========================================="
echo ""
echo "Output: $OUTPUT_DIR/eQTL_nominal_${OPTIMAL_PCS}PC_all.txt.gz"
echo ""
echo "Use this file for:"
echo "  - Colocalization analysis"
echo "  - Fine-mapping"
echo "  - Effect size estimation"
echo "  - Mendelian randomization"
echo ""
