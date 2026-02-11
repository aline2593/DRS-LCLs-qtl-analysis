# Gene Expression Quantification

This directory contains scripts for quantifying gene-level expression from Direct RNA long-reads sequencing (DRS) and preparing data for QTL analysis.

---

## Quick Start

```bash
# Complete gene expression workflow
bash 01_longread_alignment_quantification.sh
bash 02_gene_expression_filtering.sh
Rscript 03_longread_vs_shortread_comparison.R

# Then proceed to QTL mapping
cd ../03_qtl_mapping
bash 01_eqtl_permutation.sh
```

---

## Scripts Overview

**01. `01_longread_alignment_quantification.sh`**
- Aligns long-read RNA-seq to reference genome (minimap2)
- Counts reads per gene (featureCounts)
- Normalizes to RPKM (edgeR)
- Creates BED file with genomic coordinates
- **Output**: Gene expression matrix, sorted BAMs

**02. `02_gene_expression_filtering.sh`**
- Filters lowly expressed genes (50% threshold)
- Compresses and indexes for QTL analysis
- **Output**: Filtered BED file ready for QTLtools

**03. `03_longread_vs_shortread_comparison.R`**
- Compares long-read vs short-read quantification
- Correlation and concordance analysis
- Identifies technology-specific genes
- **Output**: Comparison plots and statistics

---

## Complete Workflow

### Step 1: Alignment and Quantification
```bash
bash 01_longread_alignment_quantification.sh
```

**What it does:**
1. Aligns reads with minimap2 (splice-aware, long-read mode)
2. Converts to BAM, sorts, and indexes
3. Counts reads per gene with featureCounts
4. Calculates RPKM normalization
5. Creates BED file with coordinates + expression

**Key Parameters:**
- `-ax splice` - Splice-aware alignment for long reads
- `-k14 -w5` - Minimizer parameters optimized for long reads
- `--fracOverlap 0.2` - Require 20% overlap for assignment
- `--largestOverlap` - Assign to gene with largest overlap

**Output:**
- `sorted_alignment_*.bam` - Aligned reads (one per sample)
- `sorted_alignment_*.bam.bai` - BAM indexes
- `Count_matrix_final.txt` - Raw read counts per gene
- `rpkm_expression_matrix.txt` - RPKM-normalized expression
- `gene_expression_protein_coding_lncRNA.bed` - BED format for QTL

---

### Step 2: Expression Filtering
```bash
bash 02_gene_expression_filtering.sh
```

**What it does:**
1. Loads gene expression BED file
2. Calculates proportion of zero-expressed samples per gene
3. Keeps genes expressed (RPKM > 0) in ≥50% of samples
4. Compresses with bgzip and indexes with tabix

**Why 50% threshold?**
- Balances statistical power vs false discovery
- Removes genes with sporadic/unreliable detection across samples
- Reduces multiple testing burden in QTL analysis
- Standard threshold in eQTL studies

**Output:**
- `gene_expression_filtered50.bed.gz` - Compressed, filtered expression
- `gene_expression_filtered50.bed.gz.tbi` - Index for fast access

---

### Step 3: Technology Comparison
```bash
Rscript 03_longread_vs_shortread_comparison.R
```

**What it does:**
1. Loads both long-read (nanopore) and short-read (Illumina) expression
2. Identifies overlapping genes between technologies
3. Calculates Spearman correlations (matched samples and all pairwise)
4. Compares expression levels between shared and unique genes
5. Analyzes gene length differences

**Output:**
- `correlation_heatmap_matched.pdf` - Same sample correlations
- `correlation_histogram.pdf` - Distribution of correlations
- `expression_comparison_violin.pdf` - Shared vs unique genes
- `gene_length_comparison.pdf` - Length distributions
- `technology_comparison_summary.txt` - Statistics

---

## Output File Formats

### Gene Expression BED (Unfiltered)
```
#chr  start  end    gene_id          .  strand  Sample1  Sample2  ...
chr1  11869  14409  ENSG00000223972  .  +       2.45     3.12     ...
chr1  14404  29570  ENSG00000227232  .  -       0.89     1.23     ...
```

**Columns:**
1. Chromosome
2. Gene start position
3. Gene end position
4. Gene ID (Ensembl)
5. Score (not used, always ".")
6. Strand (+/-)
7-N. RPKM expression for each sample

### Filtered BED (For QTL Analysis)
Same format, but:
- Only genes expressed in ≥50% of samples
- Compressed with bgzip
- Indexed with tabix for fast random access

---

## Dependencies

### Software
- **minimap2** (v2.17+) - Long-read aligner
- **samtools** (v1.10+) - BAM file manipulation
- **featureCounts** (subread package v2.0+) - Read counting
- **bgzip** and **tabix** (htslib) - Compression and indexing

### R Packages
- `edgeR` - RPKM normalization
- `data.table` - Fast data manipulation
- `Hmisc` - Correlation with p-values (for comparison script)
- `corrplot` - Correlation visualization (for comparison script)
- `ggplot2` - Plotting (for comparison script)

### Installation
```bash
# minimap2
conda install -c bioconda minimap2

# samtools
conda install -c bioconda samtools

# subread (featureCounts)
conda install -c bioconda subread

# htslib (bgzip, tabix)
conda install -c bioconda htslib

# R packages
R -e "install.packages('BiocManager'); BiocManager::install(c('edgeR', 'Hmisc'))"
R -e "install.packages(c('data.table', 'corrplot', 'ggplot2'))"
```

---

## Troubleshooting

**Issue: minimap2 alignment rate very low (<50%)**
→ Check that FASTQ files are long-read RNA-seq (not DNA)
→ Verify reference genome matches sample species
→ Ensure using `-ax splice` 

**Issue: featureCounts fails with "no features"**
→ Check chromosome naming: GTF and BAM must match (chr1 vs 1)
→ Use `--checkFragLength` to see if reads are assigned

**Issue: Low correlation with short-read**
→ Different normalization methods can affect correlation
→ Batch effects between sequencing runs
→ Check that correct samples are being compared

---

## Notes

### RPKM vs TPM
- This pipeline uses **RPKM** (Reads Per Kilobase Million)
- TPM (Transcripts Per Million) is more common now
- For QTL analysis, choice doesn't matter (rank-normalized anyway)
- Both account for gene length and library size

### 50% Expression Filter
- **Conservative**: Ensures reliable detection
- Alternative thresholds:
  * 30% for pilot studies or low depth
  * 70% for very high-quality datasets
- Can adjust in script 02 if needed

---

## Next Steps

After completing gene expression quantification:

1. **Proceed to QTL mapping** (recommended):
   ```bash
   cd ../03_qtl_mapping
   bash 01_eqtl_permutation.sh  # Discover eQTLs
   ```

2. **Alternative**: Analyze transcript isoforms:
   ```bash
   cd ../02_transcript_quantification
   bash 01_flair_transcript_quantification.sh
   ```

3. **Compare with short-read replication cohort**
   - Use script 03 to assess concordance
   - Important for validating nanopore quantification

---

## Citation

If using these scripts, please cite:

**Tools:**
- minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics.
- featureCounts: Liao, Y., et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics.
- edgeR: Robinson, M.D., et al. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics.

**Your study:**
https://www.researchsquare.com/article/rs-4613444/v1
