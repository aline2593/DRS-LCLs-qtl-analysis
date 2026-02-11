#!/bin/bash
################################################################################
# Transcript-Level Quantification with FLAIR
################################################################################
# Description: Quantify transcript isoforms from Direct RNA-seq data using FLAIR
#              (Full-Length Alternative Isoform analysis of RNA). FLAIR identifies
#              and quantifies high-confidence transcript isoforms by correcting,
#              collapsing, and quantifying long reads.
#
# Input:  - Sorted BAM files from minimap2 alignment
#         - Reference genome (FASTA)
#         - Gene annotation (GTF)
#         - Raw FASTQ files (for collapse step)
#
# Output: - Corrected read alignments (BED12 format)
#         - High-confidence isoform sequences (FASTA)
#         - Transcript quantification matrix (TPM values)
#
# Dependencies: FLAIR, samtools, python
#
# Reference: Tang et al. (2020) FLAIR: Full-Length Alternative Isoform
#            analysis of RNA. Nature Communications.
################################################################################

################################################################################
# SETUP: Create conda environment and install FLAIR
################################################################################

# Create FLAIR environment (one-time setup)
# mamba create -n flair -c conda-forge -c bioconda flair

# Activate FLAIR environment
conda activate flair

################################################################################
# STEP 1: INDEX BAM FILES
################################################################################
# Bam file from aligment using minimap2 from 01_longread_alignment_gene_quant.sh
# Ensure all BAM files are indexed (required for FLAIR)

echo "Indexing BAM files..."

for bam_file in /path/to/bam/*.bam; do
  if [[ -f "$bam_file" ]]; then
    samtools index "$bam_file"
    echo "  Indexed: $(basename $bam_file)"
  fi
done

echo "Indexing complete"

################################################################################
# STEP 2: CONVERT BAM TO BED12 FORMAT
################################################################################
# Convert aligned reads to BED12 format (required input for FLAIR correct)
# BED12 format includes splice junction information

echo "Converting BAM to BED12 format..."

mkdir -p 01_bed12

for bam_file in /path/to/bam/*.bam; do
  if [[ -f "$bam_file" ]]; then
    sample_name=$(basename "$bam_file" .bam)

    echo "  Processing: $sample_name"

    python /path/to/flair/bin/bam2Bed12.py \
      -i "$bam_file" \
      > 01_bed12/${sample_name}.bed12
  fi
done

echo "BED12 conversion complete"

################################################################################
# STEP 3: CORRECT READ ALIGNMENTS
################################################################################
# Correct splice junctions in reads using reference annotation
# This step:
#   - Adjusts splice sites to match known junctions in GTF
#   - Filters out reads with low-quality alignments
#   - Generates corrected BED12 files for each sample
#   -t parameter is setting the number of threads

echo "Correcting read alignments..."

mkdir -p 02_corrected

for bed_file in 01_bed12/*.bed12; do
  if [[ -f "$bed_file" ]]; then
    sample_name=$(basename "$bed_file" .bed12)

    echo "  Correcting: $sample_name"

    python /path/to/flair/bin/flair.py correct \
      -q "$bed_file" \
      -g /path/to/reference/hg19.fa \
      -f /path/to/annotation/gencode.v46lift37.annotation.gtf \
      -t 10 \
      -o 02_corrected/${sample_name}

    # Output will be: ${sample_name}_all_corrected.bed
  fi
done

echo "Correction complete"

################################################################################
# STEP 4: CONCATENATE CORRECTED READS FROM ALL SAMPLES
################################################################################
# Merge corrected reads from all samples for isoform discovery
# This creates a comprehensive set of isoforms across all samples

echo "Concatenating corrected reads..."

mkdir -p 03_concatenated

cat 02_corrected/*_all_corrected.bed > 03_concatenated/all_samples_corrected.bed

# Count total reads
total_reads=$(wc -l < 03_concatenated/all_samples_corrected.bed)
echo "Total corrected reads: $total_reads"

################################################################################
# STEP 5: SPLIT CONCATENATED FILE BY CHROMOSOME
################################################################################
# Split by chromosome to enable parallel processing of collapse step (memeory intensive)
# Processing by chromosome reduces memory requirements

echo "Splitting by chromosome..."

mkdir -p 03_concatenated/by_chromosome

# Array of chromosomes to process
chromosomes=(chr{1..22} chrX chrY chrM)

for chr in "${chromosomes[@]}"; do
  echo "  Extracting: $chr"

  # Extract reads from this chromosome (match whole word to avoid chr1 matching chr10)
  grep -w "^${chr}" 03_concatenated/all_samples_corrected.bed \
    > 03_concatenated/by_chromosome/all_samples_corrected.${chr}.bed

  # Count reads for this chromosome
  chr_reads=$(wc -l < 03_concatenated/by_chromosome/all_samples_corrected.${chr}.bed)
  echo "    Reads on $chr: $chr_reads"
done

echo "Chromosome splitting complete"

################################################################################
# STEP 6: PREPARE CONCATENATED FASTQ FOR COLLAPSE
################################################################################
# FLAIR collapse requires the original FASTQ reads to generate consensus sequences
# Concatenate all sample FASTQs into a single file

echo "Preparing concatenated FASTQ..."

# If FASTQs are in separate files per sample
zcat /path/to/fastq/*.fastq.gz | gzip > 04_collapse/all_samples_concatenated.fastq.gz

echo "FASTQ concatenation complete"

################################################################################
# STEP 7: COLLAPSE READS INTO HIGH-CONFIDENCE ISOFORMS
################################################################################
# Collapse corrected reads into unique isoform models
# This step:
#   - Groups reads with similar splice patterns
#   - Generates consensus isoform sequences
#   - Filters low-confidence isoforms
#   - Produces FASTA files of isoform sequences

echo "Collapsing reads into isoforms (by chromosome)..."

mkdir -p 04_collapse

# Process each chromosome separately
for bed_file in 03_concatenated/by_chromosome/all_samples_corrected.*.bed; do
  if [[ -f "$bed_file" ]]; then

    # Extract chromosome name from filename
    chr_name=$(basename "$bed_file" | sed 's/all_samples_corrected\.//' | sed 's/\.bed$//')

    echo "  Collapsing chromosome: $chr_name"

    python /path/to/flair/bin/flair.py collapse \
      -g /path/to/reference/hg19.fa \
      -r 04_collapse/all_samples_concatenated.fastq.gz \
      -q "$bed_file" \
      -f /path/to/annotation/gencode.v46lift37.annotation.gtf \
      -s 10 \
      --quality 1 \
      --keep_intermediate \
      --temp_dir 04_collapse/temp_${chr_name} \
      -o 04_collapse/${chr_name}

    # Output files:
    #   ${chr_name}.isoforms.fa     - Isoform sequences (FASTA)
    #   ${chr_name}.isoforms.bed    - Isoform coordinates (BED12)
    #   ${chr_name}.isoforms.gtf    - Isoform annotation (GTF)

  fi
done

# Parameters explained:
#   -s 10: Minimum number of supporting reads for an isoform
#   --quality 1: Minimum alignment quality (1 = lenient, allows more isoforms)
#   --keep_intermediate: Save intermediate files for debugging
#   --temp_dir: Directory for temporary files

echo "Collapse complete"

################################################################################
# STEP 8: PREPARE SAMPLE MANIFEST FOR QUANTIFICATION
################################################################################
# Create a manifest file mapping sample IDs to their FASTQ files
# Format: sample_id<TAB>fastq_path

echo "Creating sample manifest..."

# Example manifest content:
# sample1	/path/to/sample1.fastq.gz
# sample2	/path/to/sample2.fastq.gz
# ...

cat > 05_quantify/samples_manifest.tsv << EOF
# sample_id	fastq_path
EOF

# Generate manifest entries for all samples
for fastq_file in /path/to/fastq/*.fastq.gz; do
  sample_id=$(basename "$fastq_file" .fastq.gz)
  echo -e "${sample_id}\t${fastq_file}" >> 05_quantify/samples_manifest.tsv
done

echo "Sample manifest created: 05_quantify/samples_manifest.tsv"

################################################################################
# STEP 9: QUANTIFY ISOFORM EXPRESSION
################################################################################
# Quantify expression of discovered isoforms in each sample
# This step:
#   - Counts reads supporting each isoform per sample
#   - Normalizes to TPM (Transcripts Per Million)
#   - Outputs expression matrix

echo "Quantifying isoform expression (by chromosome)..."

mkdir -p 05_quantify

# Quantify each chromosome separately
for fasta_file in 04_collapse/*.isoforms.fa; do
  if [[ -f "$fasta_file" ]]; then

    # Extract chromosome name from filename
    chr_name=$(basename "$fasta_file" .isoforms.fa)

    echo "  Quantifying chromosome: $chr_name"

    python /path/to/flair/bin/flair.py quantify \
      -r 05_quantify/samples_manifest.tsv \
      -i "$fasta_file" \
      --tpm \
      --temp_dir 05_quantify/temp_${chr_name} \
      -o 05_quantify/${chr_name}

    # Output file:
    #   ${chr_name}.counts.tsv - Raw counts and TPM values per sample

  fi
done

# Parameters explained:
#   -r: Sample manifest (TSV format)
#   -i: Isoform sequences (FASTA from collapse step)
#   --tpm: Calculate TPM normalization
#   -o: Output prefix

echo "Quantification complete"

################################################################################
# STEP 10: MERGE QUANTIFICATION RESULTS ACROSS CHROMOSOMES
################################################################################
# Combine per-chromosome quantification into genome-wide matrix

echo "Merging quantification results..."

# Combine all chromosome-level count files
# First file includes header, subsequent files skip header
head -n 1 05_quantify/chr1.counts.tsv > 05_quantify/all_chromosomes_counts.tsv

for counts_file in 05_quantify/*.counts.tsv; do
  chr_name=$(basename "$counts_file" .counts.tsv)

  # Skip if this is chr1 (already added with header) or the merged file
  if [[ "$chr_name" != "chr1" ]] && [[ "$chr_name" != "all_chromosomes_counts" ]]; then
    tail -n +2 "$counts_file" >> 05_quantify/all_chromosomes_counts.tsv
  fi
done

# Count total isoforms
total_isoforms=$(tail -n +2 05_quantify/all_chromosomes_counts.tsv | wc -l)
echo "Total isoforms quantified: $total_isoforms"

################################################################################
# STEP 11: MERGE ISOFORM ANNOTATIONS
################################################################################
# Combine per-chromosome isoform BED and GTF files

echo "Merging isoform annotations..."

# Merge BED files
cat 04_collapse/*.isoforms.bed > 04_collapse/all_chromosomes.isoforms.bed

# Merge GTF files (skip headers)
for gtf_file in 04_collapse/*.isoforms.gtf; do
  # GTF files may have headers starting with #
  grep -v "^#" "$gtf_file"
done > 04_collapse/all_chromosomes.isoforms.gtf

# Merge FASTA files
cat 04_collapse/*.isoforms.fa > 04_collapse/all_chromosomes.isoforms.fa

echo "Annotation merging complete"

################################################################################
# SUMMARY
################################################################################

echo ""
echo "========================================="
echo "FLAIR Pipeline Complete!"
echo "========================================="
echo ""
echo "Output files:"
echo "  Corrected reads:"
echo "    02_corrected/*_all_corrected.bed"
echo ""
echo "  Isoform sequences (per chromosome):"
echo "    04_collapse/*.isoforms.fa"
echo "    04_collapse/*.isoforms.bed"
echo "    04_collapse/*.isoforms.gtf"
echo ""
echo "  Isoform sequences (genome-wide):"
echo "    04_collapse/all_chromosomes.isoforms.fa"
echo "    04_collapse/all_chromosomes.isoforms.bed"
echo "    04_collapse/all_chromosomes.isoforms.gtf"
echo ""
echo "  Expression quantification:"
echo "    05_quantify/*.counts.tsv (per chromosome)"
echo "    05_quantify/all_chromosomes_counts.tsv (genome-wide)"
echo ""
echo "Summary statistics:"
echo "  Total corrected reads: $total_reads"
echo "  Total isoforms identified: $total_isoforms"
echo ""
echo "Next steps:"
echo "  1. Filter lowly expressed isoforms"
echo "  2. Compare with known annotations"
echo "  3. Perform transcript-level QTL mapping"
echo "========================================="
