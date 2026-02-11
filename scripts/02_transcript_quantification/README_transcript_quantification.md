# Transcript Isoform Quantification from Long-Read RNA-seq

This directory contains scripts for quantifying transcript isoforms from Direct RNA long-reads sequencing data (DRS). Three different methods were tested, with **FLAIR** providing the results used in the main manuscript.

## Overview of Methods

| Method | Used in Paper | Description |
|--------|---------------|-------------|
| **FLAIR** | Main text | Full-Length Alternative Isoform analysis - correction, collapse, quantification |
| **Bambu** |Supplementary Notes| Isoform discovery and quantification using Bayesian model |
| **StringTie** | Supplementary Notes| Transcript assembly and quantification |

---

## Scripts

### 1. `01_flair_transcript_quantification.sh` => Main Method
**Purpose:** Identify and quantify high-confidence transcript isoforms using FLAIR

**Method Overview:**
FLAIR (Full-Length Alternative Isoform analysis of RNA) is a computational pipeline specifically designed for identifying isoforms from long-read RNA sequencing. It uses a three-step approach:
1. **Correct**: Adjust splice junctions to match known annotations
2. **Collapse**: Group similar reads into consensus isoform models
3. **Quantify**: Count reads supporting each isoform across samples

**Input:**
- Sorted BAM files from minimap2 alignment (from gene expression pipeline)
- Reference genome (hg19 FASTA)
- Gene annotation (Gencode v46 GTF)
- Raw FASTQ files (for consensus sequence generation)

**Output:**
- Corrected read alignments (`02_corrected/*_all_corrected.bed`)
- High-confidence isoform sequences (`04_collapse/*.isoforms.fa`)
- Isoform coordinates (`04_collapse/*.isoforms.bed`)
- Isoform annotations (`04_collapse/*.isoforms.gtf`)
- Expression matrix (`05_quantify/all_chromosomes_counts.tsv`)
  - Format: isoform_id × samples with raw counts and TPM values

**Pipeline Steps:**
1. Index BAM files
2. Convert BAM to BED12 format (splice junction format)
3. Correct splice sites using reference annotation
4. Concatenate corrected reads from all samples
5. Split by chromosome for parallel processing
6. Collapse reads into consensus isoforms (min 10 supporting reads)
7. Prepare sample manifest
8. Quantify isoform expression (TPM normalization)
9. Merge results across chromosomes

**Key Parameters:**
- `--s 10`: Minimum 10 supporting reads per isoform (balance sensitivity vs noise)
- `--quality 1`: Lenient quality threshold to capture diverse isoforms
- `--tpm`: Normalize to Transcripts Per Million

**Computational Requirements:**
- Memory: ~20-40GB per chromosome for collapse step
- Time: ~4-8 hours for 60 samples (depending on read depth)
- Storage: ~50-100GB for intermediate files

---

### 2. `02_flair_complete_analysis_with_sqanti.R` => Complete Downstream Analysis
**Purpose:** Comprehensive analysis from FLAIR quantification to publication-ready results in a single workflow

**What This Script Does:**
This consolidated script combines all downstream analysis steps:
1. Load and filter FLAIR TPM quantification
2. Run SQANTI3 for structural classification (automated)
3. Reannotate using SQANTI results
4. Filter to protein-coding/lincRNA at 50% expression
5. Generate all final publication plots

**Why Combined Workflow?**
Instead of running 3 separate scripts (initial filtering → SQANTI → reannotation), this single script handles the entire process, avoiding intermediate files and manual steps.

**Complete Workflow:**
```
Part 1: Load FLAIR quantification → Filter unknowns → TPM >= 5 filter
Part 2: Create filtered GTF → Run SQANTI3 (system call)
Part 3: Load SQANTI results → Reannotate novel isoforms
Part 4: Remove duplicate transcript IDs
Part 5: Filter to protein-coding/lincRNA only
Part 6: Create 50% filtered BED for tQTL
Part 7: Calculate expression statistics
Part 8: Generate all publication plots
```

**SQANTI3 Classification:**
SQANTI3 classifies isoforms by structural categories:
- **full-splice_match (FSM)**: All junctions match known transcript (Annotated)
- **incomplete-splice_match (ISM)**: Subset of known junctions (Novel)
- **novel_in_catalog (NIC)**: New combination of known sites (Novel)
- **novel_not_in_catalog (NNC)**: Novel splice sites (Novel)

**Reannotation Logic:**
- If initially "novel" BUT SQANTI associates with ENST → reclassify as "annotated"
- Remove duplicate ENST IDs (same transcript counted twice)
- Result: accurate classification of truly novel isoforms

**Input:**
- FLAIR TPM quantification files (`chr*.tpm.tsv`)
- Gencode annotation (GTF)
- Reference genome (FASTA)
- PolyA motif list (for SQANTI)
- Gene-level expression data

**Output:**
- `01_isoform_summary_before_sqanti.txt` - Initial classification
- `02_isoforms_reannotated_sqanti.txt` - After SQANTI reannotation
- `03_Novel_Proteincoding_lincRNA.txt` - Novel isoforms (final)
- `04_Annotated_Proteincoding_lincRNA.txt` - Annotated isoforms (final)
- `05_transcripts_filtered50_for_tQTL.bed` **Use for tQTL mapping**
- Publication plots (PDFs):
  - Transcript expression distributions
  - Sample breadth analysis
  - Gene expression vs isoform number (3 plots)
- `sqanti_output/` - Complete SQANTI3 results and QC reports

**Run Time:** ~2-4 hours total (SQANTI is the slowest step)

**Configuration Required:**
Update these paths in the script:
- `/path/to/SQANTI3/` - SQANTI3 installation
- `/path/to/annotation/gencode.v46lift37.annotation.gtf`
- `/path/to/reference/hg19.fa`
- `/path/to/polyA_motifs/mouse_and_human.polyA_motif.txt`
- `/path/to/gene/expression/rpkm_edgeR.txt`

---

### 3. `03_bambu_transcript_quantification.R` (Supplementary)
**Purpose:** Alternative isoform quantification using Bambu's Bayesian framework

**Method Overview:**
Bambu uses a Bayesian approach to discover novel isoforms while leveraging existing annotations. It's particularly good at identifying lowly expressed isoforms and provides uncertainty estimates. Paticualrly biased toward annotated isoforms.

**Input:**
- Same BAM files as FLAIR
- Reference genome and annotation

**Output:**
- Extended isoform annotations
- Expression estimates with credible intervals
- Novel vs annotated isoform classification

**Key Differences from FLAIR:**
- Bayesian statistical framework
- Provides confidence intervals
- Generally more conservative in novel isoform discovery

**Status:** Tested for comparison; results in Supplementary Notes

---

### 3. `03_stringtie_transcript_quantification.sh` (Supplementary)
**Purpose:** Transcript assembly and quantification using StringTie

**Method Overview:**
StringTie assembles transcripts from aligned reads using a network flow algorithm. Originally designed for short reads, can also work with long reads.

**Input:**
- Sorted BAM files
- Reference annotation (optional, for guided assembly)

**Output:**
- Assembled transcript GTF files
- Abundance estimates (FPKM/TPM)

**Key Differences from FLAIR:**
- De novo assembly approach
- May produce more fragmented transcripts with long reads
- Fast and memory-efficient

**Status:** Tested for comparison; results in Supplementary Notes

---

## Dependencies

### Software
- **FLAIR** (conda install)
  - Python 3.7+
  - minimap2 (for initial alignment)
  - samtools (for BAM manipulation)
- **Bambu** (R package from Bioconductor)
- **StringTie** (v2.0+)

### R packages (for downstream analysis and Bambu)
- `data.table` - Fast data reading and manipulation
- `parallel` - Parallel processing for speed
- `ggplot2` - Visualization
- `reshape2` - Data reshaping
- `gplots` - Heatmaps
- `tidyr` - Data tidying and parsing
- `bambu` - Alternative quantification method (Bioconductor)

---

## Usage

### Complete FLAIR Workflow (Recommended - All in One)

```bash
# Step 1: Run FLAIR quantification
conda activate flair
bash 01_flair_transcript_quantification.sh

# Step 2: Complete downstream analysis (filtering + SQANTI + reannotation + plots)
# This single script does everything!
Rscript 02_flair_complete_analysis_with_sqanti.R

# Step 3: Compress and index for QTL analysis
bgzip 05_transcripts_filtered50_for_tQTL.bed
tabix -p bed 05_transcripts_filtered50_for_tQTL.bed.gz
```

** Script 02 handles:
- Initial filtering
- Running SQANTI3
- Reannotation
- Final plots

No need to run multiple separate scripts.

---

### For Alternative Methods (Supplementary Notes)

```bash
# Bambu quantification
Rscript 03_bambu_transcript_quantification.R

# StringTie quantification  
bash 04_stringtie_transcript_quantification.sh
```

---

## Output File Formats

### FLAIR isoform FASTA (`.isoforms.fa`)
```
>chr1_1234_1
ATCGATCGATCG...
>chr1_1234_2
ATCGATCGATCG...
```
Each isoform has a unique ID: `chromosome_position_isoformNumber`

### FLAIR counts file (`.counts.tsv`)
```
isoform_id              sample1_counts  sample1_tpm  sample2_counts  sample2_tpm  ...
chr1_11869_14409_1      145             2.34         152             2.41         ...
chr1_14404_29570_1      89              1.45         95              1.52         ...
```

### FLAIR BED12 (`.isoforms.bed`)
Standard BED12 format with 12 columns:
1. Chromosome
2. Start position
3. End position
4. Isoform ID
5. Score
6. Strand
7. Thick start
8. Thick end
9. RGB color
10. Block count (exon count)
11. Block sizes (exon sizes)
12. Block starts (exon starts)

---

## Quality Control and Filtering

### Recommended Filters for Downstream Analysis

After quantification, filter isoforms to reduce noise:

```R
library(data.table)

# Load FLAIR counts
counts <- fread("05_quantify/all_chromosomes_counts.tsv")

# Extract TPM columns
tpm_cols <- grep("_tpm$", colnames(counts), value = TRUE)
tpm_data <- counts[, ..tpm_cols]

# Filter 1: Expression in at least 50% of samples
expressed <- apply(tpm_data, 1, function(x) sum(x > 0) / length(x) >= 0.5)

# Filter 2: Minimum average TPM threshold
min_avg_tpm <- apply(tpm_data, 1, mean) >= 1

# Combine filters
keep <- expressed & min_avg_tpm

filtered_counts <- counts[keep, ]
cat("Isoforms before filtering:", nrow(counts), "\n")
cat("Isoforms after filtering:", nrow(filtered_counts), "\n")

# Save filtered data
fwrite(filtered_counts, "05_quantify/filtered_isoforms_counts.tsv", sep = "\t")
```

---

## Method Comparison (for Supplementary Notes)

To compare the three methods:

```R
# Load outputs from all three methods
flair_counts <- fread("flair_output/all_chromosomes_counts.tsv")
bambu_counts <- fread("bambu_output/counts_transcript.txt")
stringtie_counts <- fread("stringtie_output/transcript_counts.tsv")

# Compare:
# 1. Number of isoforms detected
# 2. Overlap of detected isoforms
# 3. Correlation of expression values for shared isoforms
# 4. Novel vs annotated isoform ratios
```

---

## Troubleshooting

**Issue:** FLAIR correct fails with "no reads in BED file"
- Check that BED12 files are not empty: `wc -l 01_bed12/*.bed12`
- Verify BAM files contain aligned reads: `samtools view -c file.bam`

**Issue:** FLAIR collapse runs out of memory
- Split by chromosome (already implemented in script)
- Reduce number of input reads
- Increase memory allocation or use high-memory nodes

**Issue:** Very few isoforms in quantification output
- Check `-s` parameter (minimum support): lower from 10 to 5 for more isoforms
- Verify FASTQ file matches BAM file samples
- Check manifest file format (tab-separated, no spaces)

---

## Performance Optimization

### For large datasets (>100 samples):
1. **Process samples in batches** for correction step
2. **Use parallel processing** for collapse step:
   ```bash
   # Use GNU parallel
   parallel -j 4 'bash collapse_chr.sh {}' ::: chr{1..22} chrX chrY chrM
   ```
3. **Increase stringency** to reduce computation:
   - Increase `-s` to 15 or 20 (minimum supporting reads)
   - Increase `--quality` threshold

### For limited storage:
1. **Remove intermediate files** after each step:
   ```bash
   rm -rf 04_collapse/temp_*  # Remove temporary collapse files
   ```
2. **Compress count files**:
   ```bash
   gzip 05_quantify/*.counts.tsv
   ```

---

## Notes

### Why split by chromosome?
- **Memory efficiency**: Processing all reads at once can require >100GB RAM
- **Parallelization**: Chromosomes can be processed independently
- **Debugging**: Easier to identify and fix chromosome-specific issues

### Why concatenate all sample reads?
- **Comprehensive isoform discovery**: Each sample may have unique isoforms
- **Increased confidence**: More evidence for each isoform across samples
- **Reduced false positives**: Low-quality sample-specific isoforms filtered out

### FLAIR vs Gene-level quantification
- Gene-level: Sum of all isoforms from a gene locus
- Transcript-level: Individual isoform expression
- FLAIR is more sensitive to splicing differences
- Gene-level is more robust to technical noise

---

## Citation

If you use these scripts, please cite:

**Main method (FLAIR):**
- Tang, A.D., et al. (2020). Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. Nature Communications, 11(1), 1438.

**Alternative methods:**
- Chen, Y., et al. (2021). A systematic benchmark of Nanopore long read RNA sequencing for transcript level analysis in human cell lines. bioRxiv. (Bambu)
- Kovaka, S., et al. (2019). Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biology, 20(1), 278. (StringTie)

**Your paper:**
- https://www.researchsquare.com/article/rs-4613444/v1

---

## Next Steps

After transcript quantification:
1. Filter lowly expressed transcripts
2. Annotate novel vs known isoforms
3. Perform transcript-level QTL mapping (tQTL)

