#!/bin/bash
################################################################################
# Gene Expression Filtering and eQTL Analysis
################################################################################
# Description: Filter genes by expression threshold and perform eQTL discovery
#              using QTLtools with permutation testing to determine optimal
#              number of expression PCs to include as covariates.
#
# Input:  - Gene expression BED file (gene_expression_protein_coding_lncRNA.bed)
#         - Genotype VCF/BCF file
#
# Output: - Filtered expression BED (50% expression threshold)
#         - Expression PCs
#         - eQTL permutation results
#         - eQTL nominal pass results
#
# Dependencies: QTLtools, bgzip, tabix, R (qvalue, data.table)
################################################################################

################################################################################
# STEP 1: FILTER GENES BY EXPRESSION
################################################################################
# Keep only genes expressed (RPKM > 0) in at least 50% of samples
# This reduces multiple testing burden and focuses on reliably detected genes

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

cat("Genes after 50% expression filtering:", nrow(filtered_bed), "\n")
cat("Genes removed:", sum(!pass_filter), "\n")

# Save filtered BED file
write.table(filtered_bed,
            file = "gene_expression_filtered50.bed",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

cat("Filtered BED file created: gene_expression_filtered50.bed\n")

EOF

################################################################################
# STEP 2: COMPRESS AND INDEX BED FILE
################################################################################
# QTLtools requires bgzip-compressed and tabix-indexed BED files

bgzip gene_expression_filtered50.bed
tabix -p bed gene_expression_filtered50.bed.gz

echo "BED file compressed and indexed"

################################################################################
# STEP 3: CALCULATE EXPRESSION PRINCIPAL COMPONENTS
################################################################################
# Calculate PCs from gene expression to control for hidden confounders
# Will test different numbers of PCs (0, 3, 6, ..., 54) as covariates

module load qtltools/1.3

QTLtools pca \
  --bed gene_expression_filtered50.bed.gz \
  --scale \
  --center \
  --out gene_expression_pcs

echo "Expression PCs calculated"

################################################################################
# STEP 4: PREPARE COVARIATE FILES WITH DIFFERENT NUMBERS OF PCS
################################################################################
# Create separate covariate files with 0, 3, 6, 9, ..., 54 PCs
# First PC file has 6 header rows + 0 PCs = 6 total rows
# Each subsequent file adds 3 more PCs

# Copy the base file (with headers but no PCs)
cp gene_expression_pcs.pca merged_pcs.base

# Create files with incrementally more PCs
for num_pcs in $(seq 3 3 54); do
  total_rows=$((6 + num_pcs))  # 6 header rows + N PCs
  head -n $total_rows gene_expression_pcs.pca > merged_pcs.${num_pcs}
  echo "Created covariate file with $num_pcs PCs"
done

# Also create a 0 PC file (just headers)
head -n 6 gene_expression_pcs.pca > merged_pcs.0

echo "Covariate files created for 0, 3, 6, ..., 54 PCs"

################################################################################
# STEP 5: eQTL PERMUTATION ANALYSIS
################################################################################
# Run QTLtools in permutation mode to identify significant eGenes
# Test each PC level to determine optimal number of covariates
# Split analysis into 20 chunks for parallel processing

module load qtltools/1.3

# Run permutation pass with 0 PCs (unadjusted)
for chunk in $(seq 1 20); do
  bsub -P acc_project \
       -q premium \
       -n 1 \
       -W 20:00 \
       -R "rusage[mem=25000]" \
       -o "eQTL_perm_pc0.out" \
       -e "eQTL_perm_pc0.err" \
       "QTLtools cis \
          --vcf /path/to/genotypes.bcf \
          --bed gene_expression_filtered50.bed.gz \
          --cov merged_pcs.0 \
          --normal \
          --permute 1000 \
          --chunk $chunk 20 \
          --out gene_eQTL_perm_pc0_chunk${chunk}.txt"
done

# Run permutation pass with different numbers of PCs
for num_pcs in $(seq 3 3 54); do
  echo "Submitting permutation jobs with $num_pcs PCs"

  for chunk in $(seq 1 20); do
    bsub -P acc_project \
         -q premium \
         -n 1 \
         -W 20:00 \
         -R "rusage[mem=25000]" \
         -o "eQTL_perm_pc${num_pcs}.out" \
         -e "eQTL_perm_pc${num_pcs}.err" \
         "QTLtools cis \
            --vcf /path/to/genotypes.bcf \
            --bed gene_expression_filtered50.bed.gz \
            --cov merged_pcs.${num_pcs} \
            --normal \
            --permute 1000 \
            --chunk $chunk 20 \
            --out gene_eQTL_perm_pc${num_pcs}_chunk${chunk}.txt"
  done
done

echo "Permutation jobs submitted for PC levels: 0, 3, 6, ..., 54"

################################################################################
# STEP 6: MERGE PERMUTATION RESULTS
################################################################################
# Combine chunks for each PC level

for num_pcs in 0 $(seq 3 3 54); do
  echo "Merging chunks for $num_pcs PCs"

  # Collect all chunks for this PC level
  cat gene_eQTL_perm_pc${num_pcs}_chunk*.txt | \
    gzip -c > gene_eQTL_perm_pc${num_pcs}_merged.txt.gz

  echo "  Merged file: gene_eQTL_perm_pc${num_pcs}_merged.txt.gz"
done

################################################################################
# STEP 7: DETERMINE OPTIMAL NUMBER OF PCS
################################################################################
# Calculate FDR and pi1 statistics to find optimal PC number
# pi1 = proportion of true alternative hypotheses
# Higher pi1 and more FDR < 5% discoveries indicate better power

R --vanilla <<'EOF'

library(qvalue)
library(data.table)

# Store results for each PC level
results_summary <- data.frame()

# Create PDF with p-value histograms
pdf("eQTL_pvalue_distributions.pdf", width = 10, height = 7)

# Process each permutation result file
for (perm_file in list.files(pattern = "gene_eQTL_perm_pc.*_merged.txt.gz")) {

  cat("Processing:", perm_file, "\n")

  # Read permutation results
  # Column 20 contains the adjusted p-values from permutations
  perm_data <- fread(perm_file, data.table = FALSE)

  # Calculate q-values (FDR)
  qval_obj <- qvalue(perm_data$V20)

  # Plot p-value histogram
  hist(perm_data$V20,
       col = "grey",
       breaks = 50,
       main = paste0(perm_file, "\nπ₁ = ", round(1 - qval_obj$pi0, 3)),
       xlab = "Permutation P-value",
       ylab = "Frequency")

  # Extract number of PCs from filename
  num_pcs <- as.numeric(sub(".*_pc([0-9]+)_merged\\.txt\\.gz$", "\\1", perm_file))

  # Count significant eGenes at FDR < 5%
  num_sig <- sum(qval_obj$qvalues <= 0.05, na.rm = TRUE)

  # Store summary statistics
  results_summary <- rbind(results_summary,
                          data.frame(
                            num_pcs = num_pcs,
                            num_egenes_fdr5 = num_sig,
                            pi1 = round(1 - qval_obj$pi0, 3)
                          ))
}

dev.off()

# Sort by number of PCs
results_summary <- results_summary[order(results_summary$num_pcs), ]

# Save summary table
write.table(results_summary,
            "eQTL_PC_optimization_summary.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# Create summary plots
pdf("eQTL_PC_optimization_plots.pdf", width = 12, height = 7)

# Plot: Number of significant eGenes vs number of PCs
plot(results_summary$num_pcs,
     results_summary$num_egenes_fdr5,
     type = "b",
     pch = 19,
     col = "darkblue",
     xlab = "Number of Expression PCs",
     ylab = "Number of eGenes (FDR < 5%)",
     main = "eGene Discovery by Number of Covariate PCs")
grid()

# Plot: pi1 vs number of PCs
plot(results_summary$num_pcs,
     results_summary$pi1,
     type = "b",
     pch = 19,
     col = "darkred",
     xlab = "Number of Expression PCs",
     ylab = "π₁ (Proportion of True Alternatives)",
     main = "π₁ Statistic by Number of Covariate PCs")
grid()

dev.off()

# Print optimal number of PCs
optimal_pcs <- results_summary$num_pcs[which.max(results_summary$num_egenes_fdr5)]
cat("\n=== OPTIMAL NUMBER OF PCs ===\n")
cat("Number of PCs with maximum eGene discovery:", optimal_pcs, "\n")
cat("Number of eGenes at this level:",
    max(results_summary$num_egenes_fdr5), "\n\n")

print(results_summary)

EOF

################################################################################
# STEP 8: NOMINAL PASS WITH OPTIMAL PCS
################################################################################
# After determining optimal number of PCs, run nominal pass
# This reports all variant-gene associations (not just top per gene)
# Update the --cov parameter with your optimal PC number

OPTIMAL_PCS=6  # UPDATE THIS based on Step 7 results

module load qtltools/1.3

for chunk in $(seq 1 20); do
  bsub -P acc_project \
       -q premium \
       -n 1 \
       -W 20:00 \
       -R "rusage[mem=25000]" \
       -o "eQTL_nominal.out" \
       -e "eQTL_nominal.err" \
       "QTLtools cis \
          --vcf /path/to/genotypes.bcf \
          --bed gene_expression_filtered50.bed.gz \
          --cov merged_pcs.${OPTIMAL_PCS} \
          --normal \
          --nominal 1 \
          --chunk $chunk 20 \
          --out gene_eQTL_nominal_chunk${chunk}.txt"
done

echo "Nominal pass submitted with $OPTIMAL_PCS PCs"

# Merge nominal results
cat gene_eQTL_nominal_chunk*.txt | gzip -c > gene_eQTL_nominal_all.txt.gz

echo "Pipeline complete!"
echo "Output files:"
echo "  - Filtered expression: gene_expression_filtered50.bed.gz"
echo "  - Expression PCs: gene_expression_pcs.pca"
echo "  - Permutation results: gene_eQTL_perm_pc*_merged.txt.gz"
echo "  - PC optimization: eQTL_PC_optimization_summary.txt"
echo "  - Nominal results: gene_eQTL_nominal_all.txt.gz"
