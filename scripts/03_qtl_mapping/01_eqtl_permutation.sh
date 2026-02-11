#!/usr/bin/env bash
################################################################################
# eQTL Permutation Mapping with PC Optimization using QTLtools
# (https://qtltools.github.io/qtltools/)
################################################################################
# Description: Expression QTL (eQTL) mapping using QTLtools with systematic
#              optimization of principal component covariates.
#
# Approach:
#   1. Calculate expression PCs from gene expression
#   2. Test multiple PC levels (0, 3, 6, ..., 60) as covariates
#   3. Run permutation-based eQTL mapping (1000 permutations)
#   4. Calculate FDR and π₁ statistics
#   5. Determine optimal number of PCs
#
# Input:  - Gene expression BED (filtered, from gene quantification)
#         - Genotype BCF/VCF file
#         - Base covariate file (SEX, batch, genotype PCs)
#
# Output: - eQTL results for each PC level
#         - PC optimization plots
#         - Final significant eQTLs (FDR < 5%)
#
# Dependencies: QTLtools, htslib, R (qvalue, data.table)
################################################################################

set -euo pipefail

# Configuration
EXPRESSION_BED="/path/to/gene_expression_filtered50.bed.gz"
GENOTYPES="/path/to/genotypes.bcf"
BASE_COV="/path/to/base_covariates.txt"  # SEX, GEN, GPC1-3
OUTPUT_DIR="/path/to/eqtl/output"
HELPER_SCRIPT="/path/to/qtltools_runFDR_cis.R"

# QTLtools parameters
THREADS=1
MEMORY="25000"  # MB
CHUNK_SIZE=20
PERMUTATIONS=1000

# Job submission (adjust for your cluster)
PROJECT="acc_yourproject"
QUEUE="yout_queue"

module load qtltools/1.3
module load htslib

mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/pc_optimization"

echo "=========================================="
echo "eQTL Permutation Mapping with PC Optimization"
echo "=========================================="
echo ""
echo "Expression: $EXPRESSION_BED"
echo "Genotypes: $GENOTYPES"
echo "Output: $OUTPUT_DIR"
echo ""

################################################################################
# STEP 1: CALCULATE EXPRESSION PCs
################################################################################

echo "STEP 1: Calculating expression principal components..."
echo "----------------------------------------------------------------------"

# Check if expression file is compressed and indexed
if [[ ! -f "${EXPRESSION_BED}" ]]; then
    echo "Error: Expression BED file not found: ${EXPRESSION_BED}"
    exit 1
fi

if [[ ! -f "${EXPRESSION_BED}.tbi" ]]; then
    echo "Creating index for expression BED..."
    tabix -p bed "${EXPRESSION_BED}"
fi

# Calculate PCs from expression data
bsub -P "$PROJECT" -q "$QUEUE" -n 1 -W 20:10 \
    -o "$OUTPUT_DIR/logs/pca.out" \
    -e "$OUTPUT_DIR/logs/pca.err" \
    "QTLtools pca --bed ${EXPRESSION_BED} --scale --center --out ${OUTPUT_DIR}/gene_expression_pcs"

echo "Waiting for PCA job to complete..."
echo "Check: $OUTPUT_DIR/gene_expression_pcs.pca"
echo ""

# Wait for user to check PCA completed
read -p "Press Enter when PCA job has completed..."

################################################################################
# STEP 2: CREATE COVARIATE FILES FOR EACH PC LEVEL
################################################################################

echo ""
echo "STEP 2: Creating covariate files for PC optimization..."
echo "----------------------------------------------------------------------"

# Combine base covariates with expression PCs
cat "${BASE_COV}" > "${OUTPUT_DIR}/merged_covariates_base.txt"
cat "${OUTPUT_DIR}/gene_expression_pcs.pca" >> "${OUTPUT_DIR}/merged_covariates_base.txt"

# Remove duplicate header if present
sed -i '7d' "${OUTPUT_DIR}/merged_covariates_base.txt" 2>/dev/null || true

# Create R script to generate covariate files
cat > "$OUTPUT_DIR/create_covariate_files.R" << 'EOF'
#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

cov_full <- fread(input_file)

# Fixed covariates (always included)
fixed_rows <- c("SEX", "GEN", "GPC1", "GPC2", "GPC3")

# Expression PCs
pc_rows <- grep("^genes.*PC[0-9]+$", cov_full$Sample, value = TRUE)

cat("Total expression PCs available:", length(pc_rows), "\n")

# Create files for 0, 3, 6, ..., 60 PCs
for (i in seq(0, 60, 3)) {
  if (i == 0) {
    cov_out <- cov_full[Sample %in% fixed_rows]
  } else {
    # Take first i PCs
    pcs_to_use <- pc_rows[1:min(i, length(pc_rows))]
    cov_out <- cov_full[Sample %in% c(fixed_rows, pcs_to_use)]
  }

  out_file <- file.path(output_dir, paste0("covariates_", i, "PC.txt"))
  fwrite(cov_out, out_file, sep = "\t", quote = FALSE)

  cat("Created:", out_file, "with", nrow(cov_out), "covariates\n")
}
EOF

chmod +x "$OUTPUT_DIR/create_covariate_files.R"
Rscript "$OUTPUT_DIR/create_covariate_files.R" \
    "$OUTPUT_DIR/merged_covariates_base.txt" \
    "$OUTPUT_DIR/pc_optimization"

echo "Covariate files created in: $OUTPUT_DIR/pc_optimization/"
echo ""

################################################################################
# STEP 3: RUN PERMUTATION eQTL MAPPING FOR EACH PC LEVEL
################################################################################

echo "STEP 3: Running permutation eQTL mapping for each PC level..."
echo "----------------------------------------------------------------------"
echo "This will submit multiple jobs to the cluster"
echo ""

for pc_num in $(seq 0 3 60); do

    cov_file="${OUTPUT_DIR}/pc_optimization/covariates_${pc_num}PC.txt"

    if [[ ! -f "$cov_file" ]]; then
        echo "Warning: Covariate file not found: $cov_file"
        continue
    fi

    echo "Submitting jobs for ${pc_num} PCs..."

    # Run QTLtools in chunks for parallelization
    for chunk in $(seq 1 $CHUNK_SIZE); do

        job_name="eQTL_perm_${pc_num}PC_chunk${chunk}"
        out_file="${OUTPUT_DIR}/pc_optimization/eQTL_perm_${pc_num}PC_${chunk}_${CHUNK_SIZE}.txt"

        echo "QTLtools cis \
            --vcf ${GENOTYPES} \
            --bed ${EXPRESSION_BED} \
            --cov ${cov_file} \
            --normal \
            --permute ${PERMUTATIONS} \
            --chunk ${chunk} ${CHUNK_SIZE} \
            --out ${out_file}" | \
        bsub -P "$PROJECT" -q "$QUEUE" -n $THREADS -W 20:10 \
            -J "$job_name" \
            -o "$OUTPUT_DIR/logs/${job_name}.out" \
            -e "$OUTPUT_DIR/logs/${job_name}.err" \
            -R "rusage[mem=${MEMORY}]"

    done

    echo "  Submitted $CHUNK_SIZE chunks for ${pc_num} PCs"
done

echo ""
echo "All permutation jobs submitted!"
echo "Monitor jobs with: bjobs"
echo ""

read -p "Press Enter when all jobs have completed..."

################################################################################
# STEP 4: MERGE CHUNKS FOR EACH PC LEVEL
################################################################################

echo ""
echo "STEP 4: Merging chunks for each PC level..."
echo "----------------------------------------------------------------------"

for pc_num in $(seq 0 3 60); do

    echo "Processing ${pc_num} PCs..."

    # Collect all chunk files
    chunk_files="${OUTPUT_DIR}/pc_optimization/eQTL_perm_${pc_num}PC_*_${CHUNK_SIZE}.txt"
    merged_file="${OUTPUT_DIR}/pc_optimization/eQTL_perm_${pc_num}PC_merged.txt.gz"

    # Check if chunks exist
    if ! ls ${chunk_files} 1> /dev/null 2>&1; then
        echo "  Warning: No chunk files found for ${pc_num} PCs"
        continue
    fi

    # Merge and compress
    cat ${chunk_files} | gzip -c > "$merged_file"
    echo "  Created: $merged_file"

    # Clean up individual chunks (optional)
    # rm ${chunk_files}
done

echo "Merging complete!"
echo ""

################################################################################
# STEP 5: CALCULATE FDR AND GENERATE PC OPTIMIZATION PLOTS
################################################################################

echo "STEP 5: Calculating FDR and generating PC optimization plots..."
echo "----------------------------------------------------------------------"

cat > "$OUTPUT_DIR/pc_optimization_analysis.R" << 'EOF'
#!/usr/bin/env Rscript
library(qvalue)
library(data.table)

setwd(commandArgs(trailingOnly = TRUE)[1])

numpcs <- data.frame()

# P-value distributions
pdf(file = "eQTL_pvalue_distributions.pdf", useDingbats = FALSE, width = 12, height = 8)
par(mfrow = c(3, 4))

for (f in list.files(pattern = "eQTL_perm_.*PC_merged\\.txt\\.gz$")) {

  cat("Processing:", f, "\n")

  # Load permutation results
  perm_results <- fread(f, data.table = FALSE)

  # Extract number of PCs from filename
  pc_num <- as.numeric(sub(".*_(\\d+)PC_merged\\.txt\\.gz$", "\\1", f))

  # Calculate q-values
  q <- qvalue(perm_results$V22)

  # Plot p-value distribution
  hist(perm_results$V22,
       col = "grey",
       breaks = 50,
       main = paste0(pc_num, " PCs\nπ₁ = ", round(1 - q$pi0, 3)),
       xlab = "P-value",
       ylab = "Frequency")

  # Count significant eGenes at FDR < 5%
  n_sig <- sum(q$qvalues <= 0.05, na.rm = TRUE)

  # Store results
  numpcs <- rbind(numpcs, data.frame(
    n_pcs = pc_num,
    n_egenes_fdr5 = n_sig,
    pi1 = round(1 - q$pi0, 3)
  ))

  cat("  PCs:", pc_num, "| eGenes (FDR<5%):", n_sig, "| π₁:", round(1 - q$pi0, 3), "\n")
}

dev.off()

# Sort by PC number
numpcs <- numpcs[order(numpcs$n_pcs), ]

# PC optimization plots
pdf(file = "eQTL_PC_optimization.pdf", useDingbats = FALSE, height = 7, width = 12)

par(mfrow = c(1, 2))

# Plot 1: Number of eGenes vs PCs
plot(numpcs$n_pcs, numpcs$n_egenes_fdr5,
     type = "b",
     pch = 19,
     col = "darkblue",
     lwd = 2,
     xlab = "Number of Expression PCs",
     ylab = "Number of eGenes (FDR < 5%)",
     main = "eGene Discovery vs PC Number",
     cex.lab = 1.2,
     cex.axis = 1.1)

grid()

# Plot 2: π₁ vs PCs
plot(numpcs$n_pcs, numpcs$pi1,
     type = "b",
     pch = 19,
     col = "darkred",
     lwd = 2,
     xlab = "Number of Expression PCs",
     ylab = "π₁ (Proportion of True Positives)",
     main = "π₁ by PC Number",
     cex.lab = 1.2,
     cex.axis = 1.1)

grid()

dev.off()

# Save summary table
fwrite(numpcs, "eQTL_PC_optimization_summary.txt", sep = "\t", quote = FALSE)

cat("\n")
cat("======================================================================\n")
cat("PC OPTIMIZATION SUMMARY\n")
cat("======================================================================\n")
print(numpcs)
cat("======================================================================\n")
cat("\n")
cat("Optimal PC number: Check where eGene discovery plateaus\n")
cat("Typical range: 3-15 PCs for ~60 samples\n")
cat("\n")
EOF

chmod +x "$OUTPUT_DIR/pc_optimization_analysis.R"
Rscript "$OUTPUT_DIR/pc_optimization_analysis.R" "$OUTPUT_DIR/pc_optimization"

echo ""
echo "PC optimization analysis complete!"
echo "Check: $OUTPUT_DIR/pc_optimization/eQTL_PC_optimization.pdf"
echo ""

################################################################################
# STEP 6: EXTRACT SIGNIFICANT eQTLs WITH OPTIMAL PC NUMBER
################################################################################

echo "STEP 6: Extracting significant eQTLs..."
echo "----------------------------------------------------------------------"

# User should specify optimal PC number based on plots
read -p "Enter optimal number of PCs (e.g., 3, 6, 9): " OPTIMAL_PCS

OPTIMAL_FILE="${OUTPUT_DIR}/pc_optimization/eQTL_perm_${OPTIMAL_PCS}PC_merged.txt.gz"

if [[ ! -f "$OPTIMAL_FILE" ]]; then
    echo "Error: File not found: $OPTIMAL_FILE"
    exit 1
fi

echo "Using ${OPTIMAL_PCS} PCs as optimal"
echo "Calculating FDR and extracting significant eQTLs..."

# Run FDR calculation script
if [[ -f "$HELPER_SCRIPT" ]]; then
    Rscript "$HELPER_SCRIPT" \
        "$OPTIMAL_FILE" \
        0.05 \
        "${OUTPUT_DIR}/eQTLs_FDR5_${OPTIMAL_PCS}PCs"
else
    echo "Warning: Helper script not found: $HELPER_SCRIPT"
    echo "Please run FDR calculation manually"
fi

echo ""
echo "=========================================="
echo "eQTL Permutation Mapping Complete!"
echo "=========================================="
echo ""
echo "Results:"
echo "  PC optimization: $OUTPUT_DIR/pc_optimization/"
echo "  Significant eQTLs: ${OUTPUT_DIR}/eQTLs_FDR5_${OPTIMAL_PCS}PCs.significant.txt"
echo ""
echo "Next step: Run nominal pass with ${OPTIMAL_PCS} PCs"
echo ""
