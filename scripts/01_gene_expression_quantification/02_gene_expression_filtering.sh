#!/usr/bin/env bash
################################################################################
# Gene Expression Filtering and Preparation for QTL Analysis
################################################################################
# Description: Filter genes by expression threshold and prepare BED file
#              for downstream QTL mapping.
#
# Input:  - Gene expression BED file (gene_expression_protein_coding_lncRNA.bed)
#
# Output: - Filtered expression BED (50% expression threshold)
#         - Compressed and indexed for QTL analysis
#
# Dependencies: R (data.table), bgzip, tabix
################################################################################

set -euo pipefail

echo "=========================================="
echo "Gene Expression Filtering"
echo "=========================================="
echo ""

################################################################################
# STEP 1: FILTER GENES BY EXPRESSION
################################################################################
# Keep only genes expressed (RPKM > 0) in at least 50% of samples
# This reduces multiple testing burden and focuses on reliably detected genes

echo "STEP 1: Filtering genes by expression threshold..."
echo "----------------------------------------------------------------------"

R --vanilla <<'EOF'

library(data.table)

# Read expression BED file
# Structure: #chr | start | end | gene_id | . | strand | sample1 | sample2 | ...
bed_data <- fread("gene_expression_protein_coding_lncRNA.bed", data.table = FALSE)

cat("Total genes before filtering:", nrow(bed_data), "\n")

# Extract sample IDs (column 7 onwards)
sample_cols <- 7:ncol(bed_data)
sample_ids <- colnames(bed_data)[sample_cols]
cat("Number of samples:", length(sample_ids), "\n")

# Extract annotation columns (first 6 columns)
annotation <- bed_data[, 1:6]

# Extract expression values
expression <- bed_data[, sample_cols]

# Calculate proportion of zeros per gene
# For each gene (row), count how many samples have zero expression
zero_proportion <- apply(expression, 1, function(x) sum(x == 0) / length(x))

cat("Range of zero proportions:", range(zero_proportion), "\n")

# Filter: keep genes with zeros in ≤50% of samples (expressed in ≥50%)
pass_filter <- zero_proportion <= 0.5

# Create filtered dataset
filtered_bed <- cbind(annotation[pass_filter, ], expression[pass_filter, ])

cat("\nFiltering Results:\n")
cat("  Genes before filtering:", nrow(bed_data), "\n")
cat("  Genes after 50% filter:", nrow(filtered_bed), "\n")
cat("  Genes removed:", sum(!pass_filter), "\n")
cat("  Percentage kept:", round(100 * nrow(filtered_bed) / nrow(bed_data), 1), "%\n")

# Save filtered BED file
write.table(filtered_bed,
            file = "gene_expression_filtered50.bed",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

cat("\nFiltered BED file created: gene_expression_filtered50.bed\n")

EOF

echo ""

################################################################################
# STEP 2: COMPRESS AND INDEX BED FILE
################################################################################
# QTLtools requires bgzip-compressed and tabix-indexed BED files

echo "STEP 2: Compressing and indexing BED file..."
echo "----------------------------------------------------------------------"

# Check if tools are available
if ! command -v bgzip &> /dev/null; then
    echo "Error: bgzip not found. Please install htslib"
    exit 1
fi

if ! command -v tabix &> /dev/null; then
    echo "Error: tabix not found. Please install htslib"
    exit 1
fi

# Compress
bgzip -f gene_expression_filtered50.bed
echo "  Compressed: gene_expression_filtered50.bed.gz"

# Index
tabix -p bed gene_expression_filtered50.bed.gz
echo "  Indexed: gene_expression_filtered50.bed.gz.tbi"

echo ""

################################################################################
# FINAL SUMMARY
################################################################################

echo "=========================================="
echo "Gene Expression Filtering Complete!"
echo "=========================================="
echo ""
echo "Output Files:"
echo "  gene_expression_filtered50.bed.gz     - Filtered expression (compressed)"
echo "  gene_expression_filtered50.bed.gz.tbi - Index file"
echo ""
echo "Next Steps:"
echo "  1. Run eQTL analysis:"
echo "     cd ../03_qtl_mapping"
echo "     bash 01_eqtl_permutation.sh"
echo ""
echo "  2. Or compare with short-read data:"
echo "     Rscript 03_longread_vs_shortread_comparison.R"
echo ""
