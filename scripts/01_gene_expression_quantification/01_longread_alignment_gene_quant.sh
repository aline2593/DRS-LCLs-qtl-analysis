#!/bin/bash
################################################################################
# Long-read RNA-seq Alignment and Gene Expression Quantification
################################################################################
# Description: This script performs alignment of Direct RNA-seq (nanopore) reads
#              to the human reference genome, followed by gene-level quantification
#              using featureCounts and RPKM normalization.
#
# Input:  - FASTQ files (*_all.fastq.gz)
#         - Reference genome (hg19)
#         - Gencode v46 annotation (GTF)
#
# Output: - Sorted, indexed BAM files
#         - Gene count matrix
#         - RPKM-normalized expression matrix (BED format)
#
# Dependencies: minimap2, samtools, featureCounts, R (edgeR, data.table)
################################################################################

################################################################################
# STEP 1: ALIGNMENT USING MINIMAP2
################################################################################
# Align nanopore Direct RNA-seq reads to reference genome using minimap2
# Parameters:
#   -ax splice: use spliced alignment mode for RNA-seq
#   -uf: use forward strand only (for Direct RNA-seq)
#   -k14: use k-mer size of 14 (recommended for noisy long reads)

for f in *.fastq.gz
do
  # Extract sample ID from filename
  sample_id=$(echo $f | sed 's/_all.fastq.gz//')
  echo "Processing sample: $sample_id"

  fastq_file=/path/to/fastq/${sample_id}_all.fastq.gz

  sbatch --partition=shared-bigmem \
         --mem=40G \
         --time=12:00:00 \
         --wrap="minimap2 -ax splice -uf -k14 \
                /path/to/reference/hg19.fa \
                $fastq_file > alignment_${sample_id}.sam"
done

################################################################################
# STEP 2: CONVERT SAM TO BAM
################################################################################
# Convert SAM files to compressed BAM format to save storage space

for f in *.sam
do
  sample_id=$(echo $f | sed 's/alignment_//' | sed 's/.sam//')
  echo "Converting to BAM: $sample_id"

  sam_file=/path/to/sam/alignment_${sample_id}.sam

  sbatch --partition=shared-cpu \
         --mem=10G \
         --time=12:00:00 \
         --wrap="samtools view -S -b $sam_file > alignment_${sample_id}.bam"
done

################################################################################
# STEP 3: SORT BAM FILES BY COORDINATE
################################################################################
# Sort BAM files by genomic coordinates (required for indexing and quantification)

for f in *.sam
do
  sample_id=$(echo $f | sed 's/alignment_//' | sed 's/.sam//')
  echo "Sorting BAM: $sample_id"

  sam_file=/path/to/sam/alignment_${sample_id}.sam

  samtools sort -o sorted_alignment_${sample_id}.bam $sam_file
done

################################################################################
# STEP 4: INDEX BAM FILES
################################################################################
# Create BAM index files (.bai) for rapid access

for f in sorted_alignment_*.bam
do
  echo "Indexing: $f"
  samtools index $f
done

################################################################################
# STEP 5: GENE-LEVEL QUANTIFICATION WITH FEATURECOUNTS
################################################################################
# Count reads mapping to each gene using Gencode v46 annotation
# Parameters:
#   -a: annotation file (GTF)
#   -t exon: count at exon level
#   -g gene_id: summarize by gene_id
#   --minOverlap 10: minimum 10bp overlap required
#   -M: count multi-mapping reads
#   --fraction: fractional counting for multi-mapping reads
#   -L: count long reads (important for nanopore data)

for f in /path/to/bam/sorted_alignment_*.bam
do
  # Extract sample ID
  sample_id=$(basename "$f" | sed 's/sorted_alignment_//' | sed 's/.bam//')
  echo "Quantifying genes for: $sample_id"

  featureCounts \
    -a /path/to/annotation/gencode.v46lift37.annotation.gtf \
    -t exon \
    -g gene_id \
    --minOverlap 10 \
    -M \
    --fraction \
    -L \
    -o /path/to/output/${sample_id}_gene_counts.txt \
    $f
done

################################################################################
# STEP 6: CREATE COUNT MATRIX AND CALCULATE RPKM
################################################################################
# Combine individual count files into a single matrix and normalize to RPKM
# (Reads Per Kilobase per Million mapped reads)

# Extract gene lengths from first sample (column 6 of featureCounts output)
cut -f1,6 /path/to/output/sample1_gene_counts.txt > Gene_length.txt

# Extract gene IDs as row names
cut -f1 /path/to/output/sample1_gene_counts.txt | sed 1d > gene_ids.txt
sed -i "1i\Geneid" gene_ids.txt

# Extract count column (column 7) from each sample file
# Skip first two rows (header and comment line from featureCounts)
for f in /path/to/output/*_gene_counts.txt; do
  sample_name=$(basename "$f" | sed 's/_gene_counts.txt//')
  (echo "$sample_name"; cut -f7 "$f" | sed '1,2d') > "${f}.counts"
done

# Combine all count columns into a single matrix
# Structure: gene_id | sample1 | sample2 | ... | sampleN | gene_length
paste -d"\t" gene_ids.txt *_gene_counts.txt.counts > Count_matrix.txt

# Add gene length as last column
cut -f2 Gene_length.txt > Gene_length_only.txt
paste -d"\t" Count_matrix.txt Gene_length_only.txt > Count_matrix_with_length.txt

# Remove last row if it contains summary statistics from featureCounts
sed '$d' Count_matrix_with_length.txt > Count_matrix_final.txt

################################################################################
# STEP 7: RPKM NORMALIZATION IN R
################################################################################

R --vanilla <<'EOF'

# Load required library
library(edgeR)
library(data.table)

# Set working directory
workDir <- "/path/to/output"
setwd(workDir)

# Read count matrix
# Structure: Geneid | Sample1 | Sample2 | ... | SampleN | Gene_Length
counts_data <- read.table("Count_matrix_final.txt", header = TRUE)

# Read gene lengths
gene_lengths <- read.table("Gene_length.txt", header = TRUE)

# Separate counts and gene IDs
gene_ids <- counts_data[, 1]
count_matrix <- counts_data[, -c(1, ncol(counts_data))]  # Remove gene ID and length columns
lengths_vector <- gene_lengths[, 2]

# Calculate RPKM using edgeR
# RPKM = (reads mapped to gene * 10^9) / (total reads * gene length in bp)
rpkm_matrix <- rpkm(count_matrix, lengths_vector)

# Add gene IDs as row names
rownames(rpkm_matrix) <- gene_ids

# Add gene IDs as first column for output
rpkm_output <- cbind(Geneid = gene_ids, rpkm_matrix)

# Save RPKM matrix
write.table(rpkm_output,
            file = "rpkm_expression_matrix.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

print("RPKM calculation complete!")

EOF

################################################################################
# STEP 8: CREATE BED FILE AND FILTER FOR PROTEIN-CODING AND LINCRNA
################################################################################
# Convert RPKM matrix to BED format with genomic coordinates
# BED format: #chr | start | end | gene_id | . | strand | sample1 | sample2 | ...

R --vanilla <<'EOF'

library(data.table)

# Read Gencode annotation to extract gene coordinates
gtf_file <- "/path/to/annotation/gencode.v46lift37.annotation.gtf"
gtf_data <- fread(gtf_file, header = FALSE, sep = "\t")

# Filter for gene features only
genes_data <- gtf_data[V3 == "gene"]

# Extract relevant fields from GTF
# Structure: chr | feature | start | end | strand | gene_id | gene_type
genes_annotation <- genes_data[, .(
  chr = V1,
  start = V4,
  end = V5,
  strand = V7,
  gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
  gene_type = sub('.*gene_type "([^"]+)".*', '\\1', V9)
)]

# Filter for protein-coding genes and lincRNAs only
genes_filtered <- genes_annotation[gene_type %in% c("protein_coding", "lincRNA")]

# Add placeholder column for BED format
genes_filtered$gid <- "."

# Read RPKM expression data
rpkm_data <- fread("rpkm_expression_matrix.txt", header = TRUE)

# Create BED structure: chr | start | end | gene_id | . | strand
bed_coords <- genes_filtered[, .(chr, start, end, gene_id, gid, strand)]

# Merge coordinates with expression data
bed_with_expression <- merge(bed_coords, rpkm_data, by.x = "gene_id", by.y = "Geneid")

# Reorder columns: #chr | start | end | gene_id | . | strand | samples...
bed_final <- bed_with_expression[, c(2:4, 1, 5, 6:ncol(bed_with_expression))]

# Sort by chromosome and position
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
bed_final$chr <- factor(bed_final$chr, levels = chr_order)
bed_sorted <- bed_final[order(chr, start, end)]

# Rename first column for BED format
setnames(bed_sorted, "chr", "#chr")

# Save BED file
fwrite(bed_sorted,
       "gene_expression_protein_coding_lncRNA.bed",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
       col.names = TRUE)

print("BED file created successfully!")
print(paste("Total genes:", nrow(bed_sorted)))

EOF

echo "Pipeline complete!"
echo "Output files:"
echo "  - Sorted BAM files: sorted_alignment_*.bam"
echo "  - Gene counts: *_gene_counts.txt"
echo "  - RPKM matrix: rpkm_expression_matrix.txt"
echo "  - BED file: gene_expression_protein_coding_lncRNA.bed"
