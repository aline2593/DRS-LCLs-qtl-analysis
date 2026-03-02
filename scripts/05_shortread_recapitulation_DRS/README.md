# QTL Replication Analyses

Cross-platform replication of QTL discoveries between short-read (Illumina) and long-read (Nanopore) RNA-seq in LCLs.

## Contents

- **eQTL_replication/** - Gene-level eQTL replication
- **trQTL_replication/** - Transcript-level trQTL replication
- **sQTL_replication/** - GTEx sQTL replication in nanopore trQTLs

---

## Quick Start

### 1. eQTL Replication
```bash
cd eQTL_replication/
Rscript eQTL_replication_analysis.R
```

**Input**: Short-read eQTLs → Long-read eQTLs (gene-level)

**Output**: `results/02_longread_eQTLs_replicated_FDR5.txt`

---

### 2. trQTL Replication
```bash
cd trQTL_replication/
Rscript trQTL_replication_analysis.R
```

**Input**: Short-read eQTLs → Long-read trQTLs (transcript-level)

**Output**:
- `results/03_longread_trQTLs_replicated_by_gene_FDR5.txt` (conservative)
- `results/04_longread_trQTLs_replicated_by_transcript_FDR5.txt` (liberal)

---

### 3. sQTL Replication
```bash
cd sQTL_replication/
bash filter_nanopore_trqtls.sh  # Preprocess large file
Rscript sQTL_replication_analysis.R
```

**Input**: GTEx sQTLs → Long-read trQTLs

**Output**:
- `results/05_GTEx_sQTL_replicated_in_nanopore_FDR5.txt`
- `plots/01_SQANTI3_sQTL_trQTL_structural_categories.pdf`

---

## Required Files

Each analysis directory needs its input files. See individual scripts for details.

**Common requirements**:
- Short-read eQTL results: `shortread_eQTLs_significant.txt`
- Long-read nominal results: Pre-filtered to discovery gene-SNP pairs
- SQANTI3 classifications (sQTL only): `merged_transcripts_flair_collapse_classification.txt.gz`

---

## Dependencies
```r
install.packages(c("data.table", "dplyr", "ggplot2"))
BiocManager::install(c("qvalue", "VennDiagram"))
```

---

## Key Metrics

- **π₁** (1 - π₀): Proportion of true signals
- **Replication rate**: % of discoveries replicated at FDR < 5%
- **Structural categories** (sQTL): Novel vs known isoform enrichment

---

## Analysis Comparison

| Analysis | Level | Discovery | Replication | Key Output |
|----------|-------|-----------|-------------|------------|
| eQTL | Gene | Illumina eQTLs | Nanopore eQTLs | Platform concordance |
| trQTL | Transcript | Illumina eQTLs | Nanopore trQTLs | 2 approaches (min-p vs all) |
| sQTL | Intron | GTEx sQTLs | Nanopore trQTLs | SQANTI3 structure analysis |
