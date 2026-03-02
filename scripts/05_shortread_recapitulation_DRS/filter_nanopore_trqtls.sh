#!/bin/bash
################################################################################
# Filter nanopore trQTL nominal results to GTEx sQTL gene-SNP pairs
################################################################################

# Input files (user should provide path to their nominal results)
GENE_SNP_PAIRS="results/04_GTEx_gene_SNP_pairs_for_filtering.txt"
NANOPORE_NOMINALS="nominals_norm_all_combined.noheader.txt.gz"

# Output
OUTPUT="trQTL_nominals_matching_sQTL.txt"

# Check if input files exist
if [ ! -f "$GENE_SNP_PAIRS" ]; then
    echo "ERROR: $GENE_SNP_PAIRS not found!"
    echo "Please run the R script first to generate gene-SNP pairs."
    exit 1
fi

if [ ! -f "$NANOPORE_NOMINALS" ]; then
    echo "ERROR: $NANOPORE_NOMINALS not found!"
    echo "Please provide the nanopore nominal results file."
    echo "Expected format: QTLtools nominal output (gzipped)"
    exit 1
fi

echo "Filtering nanopore trQTLs to GTEx sQTL gene-SNP pairs..."
echo "This may take several minutes..."

# Filter: match on gene (col 1, no version) and SNP rsID (col 10)
zcat "$NANOPORE_NOMINALS" | \
  awk 'BEGIN{FS=OFS=" "}
       NR==FNR {
         key[$1"\t"$2]=1  # Load gene-SNP pairs
         next
       }
       {
         gene=$1
         sub(/\..*/, "", gene)  # Remove version from gene ID
         snp=$10
         if ((gene"\t"snp) in key) print $0
       }' \
  "$GENE_SNP_PAIRS" - \
  > "$OUTPUT"

# Compress
gzip -f "$OUTPUT"

echo "✓ Complete! Output: ${OUTPUT}.gz"
echo "Lines: $(zcat ${OUTPUT}.gz | wc -l)"
