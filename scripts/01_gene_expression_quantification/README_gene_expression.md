# Gene Expression Quantification from Long-Read RNA-seq

This directory contains scripts for quantifying gene expression from Direct RNA sequencing (nanopore) data and performing expression QTL (eQTL) mapping.

## Scripts Overview

### 1. `01_longread_alignment_gene_quant.sh`
**Purpose:** Align long-read RNA-seq data and quantify gene-level expression

**Input:**
- Raw FASTQ files from Direct RNA sequencing (`*_all.fastq.gz`)
- Human reference genome (hg19)
- Gencode v46 annotation GTF file

**Output:**
- Sorted and indexed BAM files (`sorted_alignment_*.bam`)
- Gene count matrix (`Count_matrix_final.txt`)
- RPKM-normalized expression matrix (`rpkm_expression_matrix.txt`)
- BED file with expression values (`gene_expression_protein_coding_lncRNA.bed`)

**Pipeline steps:**
1. Align reads to genome using minimap2 (splice-aware, forward strand)
2. Convert SAM to BAM format
3. Sort BAM files by genomic coordinates
4. Index BAM files
5. Count reads per gene using featureCounts
6. Calculate RPKM normalization
7. Create BED file with genomic coordinates for QTL analysis

**Key parameters:**
- Only protein-coding genes and lincRNAs are retained
- Fractional counting for multi-mapping reads
- Long-read mode enabled in featureCounts

**Run time:** ~12-24 hours for 60 samples (depending on read depth)

---

### 2. `02_gene_expression_filtering_eQTL.sh`
**Purpose:** Filter lowly expressed genes and perform eQTL mapping with covariate optimization

**Input:**
- Gene expression BED file from script 01
- Genotype data (BCF/VCF format)

**Output:**
- Filtered expression BED (`gene_expression_filtered50.bed.gz`)
- Expression principal components (`gene_expression_pcs.pca`)
- eQTL permutation results for each PC level
- PC optimization summary (`eQTL_PC_optimization_summary.txt`)
- Nominal eQTL results (`gene_eQTL_nominal_all.txt.gz`)

**Pipeline steps:**
1. Filter genes: keep only those expressed (RPKM > 0) in ≥50% of samples
2. Compress and index BED file for QTLtools
3. Calculate expression PCs to control for hidden variation
4. Test different numbers of PCs (0, 3, 6, ..., 54) as covariates
5. Run eQTL permutation analysis for each PC level
6. Calculate FDR and π₁ statistics to determine optimal PC number
7. Run nominal pass with optimal number of PCs

**Key features:**
- Systematic optimization of covariate number
- Parallel processing using 20 chunks
- FDR control at 5% level
- π₁ statistic to assess proportion of true eQTL signals

**Run time:** ~2-4 hours per PC level (permutation mode)

---

### 3. `03_longread_shortread_comparison.R`
**Purpose:** Compare gene expression quantified from long-read vs short-read RNA-seq in the same samples

**Input:**
- Filtered nanopore expression BED (from script 02)
- Filtered Illumina/short-read expression BED
- Gene length annotations

**Output:**
- Correlation heatmaps and histograms
- Gene expression level comparisons (violin plots)
- Gene length comparisons
- Summary statistics table

**Analysis steps:**
1. Load and standardize gene IDs between technologies
2. Identify overlapping genes detected in both technologies
3. Calculate Spearman correlations:
   - Matched samples (same sample, different technology)
   - All pairwise correlations
   - Within-technology correlations
4. Compare expression levels:
   - Genes detected in both technologies vs technology-specific
   - Statistical testing (Wilcoxon rank-sum test)
5. Compare gene lengths between groups
6. Generate summary statistics and visualizations

**Key findings:**
- Quantifies concordance between long-read and short-read quantification
- Identifies genes uniquely detected by each technology
- Reveals expression and length characteristics of technology-specific genes

**Run time:** ~10-30 minutes (depending on number of genes)

---

## Dependencies

### Software
- **minimap2** (v2.17 or later) - Long-read aligner
- **samtools** (v1.10 or later) - BAM file manipulation
- **featureCounts** (subread package) - Read counting
- **QTLtools** (v1.3) - QTL mapping
- **bgzip** and **tabix** (htslib) - File compression and indexing

### R packages
- `edgeR` - RPKM normalization
- `data.table` - Fast data manipulation
- `qvalue` - FDR estimation for eQTL analysis
- `Hmisc` - Correlation analysis with p-values
- `corrplot` - Correlation matrix visualization
- `ggplot2` - Data visualization
- `ggsignif` - Statistical annotations on plots
- `tidyr` - Data manipulation

### System requirements
- 40GB RAM for alignment
- 25GB RAM for QTL mapping
- Access to HPC cluster with job scheduler (LSF/SLURM)

---

## Usage

### Before running:
1. Update all file paths in the scripts:
   - `/path/to/fastq/` → your FASTQ directory
   - `/path/to/reference/` → your reference genome location
   - `/path/to/annotation/` → your GTF file location
   - `/path/to/genotypes.bcf` → your genotype file
   - `/path/to/output/` → your output directory

2. Adjust job submission parameters for your cluster:
   - Change `bsub` to `sbatch` if using SLURM
   - Modify partition/queue names
   - Adjust memory and time limits

### Running the pipeline:

```bash
# Step 1: Align and quantify gene expression
bash 01_longread_alignment_gene_quant.sh

# Step 2: Filter and run eQTL analysis
bash 02_gene_expression_filtering_eQTL.sh

# After Step 2, examine PC optimization results:
cat eQTL_PC_optimization_summary.txt
# Update OPTIMAL_PCS variable in script based on results
# Re-run nominal pass section with optimal PC number

# Step 3: Compare with short-read data
Rscript 03_longread_shortread_comparison.R
```

---

## Output File Formats

### Gene expression BED file
```
#chr    start    end    gene_id         .    strand    Sample1    Sample2    ...
chr1    11869    14409  ENSG00000223972  .    +         2.45       3.12       ...
chr1    14404    29570  ENSG00000227232  .    -         0.89       1.23       ...
```

### eQTL permutation results (from QTLtools)
Key columns:
- Column 1: Phenotype ID (gene)
- Column 2: Chromosome
- Column 3-4: Phenotype start-end
- Column 8: Number of variants tested in cis
- Column 12: SNP ID of top variant
- Column 19: Beta effect size
- Column 20: Adjusted p-value (from permutations)

---

## Notes

- **Expression filtering:** The 50% threshold balances power and false discovery rate. Lowly expressed genes are removed to reduce multiple testing burden.

- **PC optimization:** More PCs control for more confounders but may also remove true genetic effects. The optimal number maximizes eGene discovery while controlling for technical variation.

- **Cis-window:** QTLtools default is 1Mb window around TSS. Adjust with `--window` parameter if needed.

- **Multi-mapping reads:** Long reads often map to multiple locations (e.g., across gene families). Fractional counting distributes these reads proportionally.

---

## Expected Results

From ~60 LCL samples with Direct RNA-seq:
- ~20,000 protein-coding genes + lincRNAs before filtering
- ~15,000-16,000 genes after 50% expression filter
- ~3,000-5,000 eGenes at FDR < 5% (depending on sample size and PC number)
- Optimal PC number typically between 3-15 for this sample size

**Comparison with short-reads:**
- ~10,000-11,000 genes overlapping between technologies
- Matched sample correlation (median): 0.85-0.92
- Technology-specific genes tend to be:
  - Lower expressed
  - Shorter in length (for nanopore-specific genes)
  - Longer in length (for Illumina-specific genes, may include poorly spliced reads)

---

## Troubleshooting

**Issue:** featureCounts fails with "no features"
- Check that GTF file chromosome names match BAM file (e.g., "chr1" vs "1")

**Issue:** QTLtools fails to find variants
- Ensure genotype file is properly indexed (`.csi` or `.tbi` file exists)
- Check chromosome naming consistency between BED and VCF/BCF

**Issue:** Low eGene discovery
- Check genotype quality and MAF filters
- Verify sample IDs match between expression and genotype files
- Consider adjusting cis-window size

---

## Citation

If you use these scripts, please cite:
- [Your paper citation here]

And the tools used:
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics.
- Liao, Y., et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics.
- Delaneau, O., et al. (2017). A complete tool set for molecular QTL discovery and analysis. Nature Communications.
