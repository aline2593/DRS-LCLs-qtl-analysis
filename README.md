[README_CLEAN_WORKFLOW.md](https://github.com/user-attachments/files/25249652/README_CLEAN_WORKFLOW.md)
# Direct RNA Sequencing QTL Analysis Pipeline

[![DOI](https://img.shields.io/badge/DOI-10.xxxx%2Fxxxxxx-blue)](https://doi.org/10.xxxx/xxxxxx)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)](https://www.linux.org/)

> A comprehensive computational pipeline for expression and transcript-level QTL discovery from Oxford Nanopore Direct RNA sequencing in lymphoblastoid cell lines

---

## Overview

This repository provides the complete analysis pipeline for discovering:
- **Expression QTLs (eQTLs)** - genetic variants affecting gene expression
- **Transcript QTLs (tQTLs)** - genetic variants affecting isoform usage  
- **m⁶A modification QTLs (m⁶AQTLs)** - genetic variants affecting RNA modifications

**Dataset:** 60 lymphoblastoid cell lines with Direct RNA sequencing and whole genome sequencing

---

## Analysis Workflow

### Stage 1: Gene Expression Quantification

```
Raw FASTQ files
      ↓
[Align to genome] ────────────── minimap2
      ↓
[Count per gene] ────────────── featureCounts  
      ↓
[Normalize to RPKM] ──────────── edgeR
      ↓
[Filter: 50% threshold]
      ↓
Output: ~15,000 genes
```

### Stage 2: Transcript Isoform Quantification

```
Aligned BAM files
      ↓
[Correct splice sites] ────────── FLAIR correct
      ↓
[Collapse to isoforms] ────────── FLAIR collapse
      ↓
[Quantify expression] ─────────── FLAIR quantify
      ↓
[Validate structure] ──────────── SQANTI3
      ↓
[Reannotate novel isoforms]
      ↓
[Filter: protein-coding/lincRNA + 50%]
      ↓
Output: ~40,000 transcripts across ~10,000 genes
```

### Stage 3: QTL Discovery

```
Expression data + Genotypes
      ↓
[Calculate expression PCs]
      ↓
[Test PC levels: 0, 3, 6, ..., 60] ──── Optimize covariates
      ↓
      ├─────────────────┬─────────────────┐
      ↓                 ↓                 ↓
[eQTL Permutation] [tQTL Permutation] [m⁶AQTL Permutation]
   1000 perms         1000 perms         1000 perms
      ↓                 ↓                 ↓
[Calculate FDR]    [Calculate FDR]    [Calculate FDR]
      ↓                 ↓                 ↓
[Nominal pass]     [Nominal pass]     [Nominal pass]
  All variants       All variants       All variants
      ↓                 ↓                 ↓
~5,000 eGenes      ~150 eTranscripts  ~800 m⁶A sites
```

### Stage 4: Replication and Validation

```
Discovery Results
      ↓
[Test in short-read cohort]
      ↓
[Compare technologies]
      ↓
[Functional annotation]
      ↓
Final validated QTLs
```

---

## Repository Structure

```
DRS-LCLs-qtl-analysis/
│
├── scripts/
│   │
│   ├── 01_gene_expression_quantification/
│   │   ├── 01_longread_alignment_quantification.sh
│   │   ├── 02_gene_expression_filtering.sh
│   │   ├── 03_longread_vs_shortread_comparison.R
│   │   └── README.md
│   │
│   ├── 02_transcript_quantification/
│   │   ├── 01_flair_transcript_quantification.sh
│   │   ├── 02_flair_complete_analysis_with_sqanti.R
│   │   ├── 03_transcript_dominance_assessment.R
│   │   ├── 04_qtl_example_plots.R
│   │   ├── 05_bambu_transcript_quantification.R
│   │   ├── 06_stringtie_transcript_quantification.sh
│   │   └── README.md
│   │
│   ├── 03_qtl_mapping/
│   │   ├── 01_eqtl_permutation.sh
│   │   ├── 02_eqtl_nominal.sh
│   │   ├── 03_tqtl_permutation.sh
│   │   ├── 04_tqtl_nominal.sh
│   │   └── README.md
│   │
│   ├── 04_m6a_modification_analysis/
│   │   ├── 01_m6a_detection.sh
│   │   ├── 02_m6a_quantification.R
│   │   ├── 03_m6a_qtl_mapping.sh
│   │   └── README.md
│   │
│   └── 05_replication_analysis/
│       ├── 01_shortread_replication.sh
│       ├── 02_compare_qtl_discoveries.R
│       └── README.md
│
├── environment/
│   ├── requirements.txt
│   └── environment.yml
│
├── README.md
└── LICENSE
```

---

## Quick Start

### Prerequisites

```bash
# Required software
minimap2 (v2.17+)
samtools (v1.10+)
featureCounts (subread v2.0+)
FLAIR (latest)
SQANTI3 (v3.0+)
QTLtools (v1.3+)
R (v4.0+) with packages: edgeR, data.table, qvalue, ggplot2
```

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/DRS-LCLs-qtl-analysis.git
cd DRS-LCLs-qtl-analysis

# Create conda environment
conda env create -f environment/environment.yml
conda activate drs-qtl
```

### Run Complete Pipeline

```bash
# Stage 1: Gene expression
cd scripts/01_gene_expression_quantification
bash 01_longread_alignment_quantification.sh
bash 02_gene_expression_filtering.sh

# Stage 2: Transcript isoforms
cd ../02_transcript_quantification
bash 01_flair_transcript_quantification.sh
Rscript 02_flair_complete_analysis_with_sqanti.R

# Stage 3: QTL mapping
cd ../03_qtl_mapping
bash 01_eqtl_permutation.sh    # Gene-level QTLs
bash 03_tqtl_permutation.sh    # Transcript-level QTLs

# Stage 4: m⁶A modifications
cd ../04_m6a_modification_analysis
bash 01_m6a_detection.sh
```

---

## Results Summary

### QTL Discovery

| Analysis Type | Features Tested | Discoveries (FDR < 5%) | Effect Size Range |
|---------------|----------------|------------------------|-------------------|
| **Gene eQTL** | 15,128 genes | 4,892 eGenes | β: -2.5 to 2.8 |
| **Transcript tQTL** | 40,584 transcripts | 187 eTranscripts | β: -1.8 to 2.1 |
| **m⁶A Modification** | 10,234 sites | 823 sites | β: -0.6 to 0.8 |

### Performance Metrics

| Step | Run Time | Memory | Storage |
|------|----------|--------|---------|
| Gene quantification | 2-4 hours/sample | 40 GB | 10 GB/sample |
| Transcript quantification | 4-8 hours total | 25 GB | 50 GB |
| eQTL permutation | 2-4 hours/PC level | 25 GB | 5 GB |
| tQTL permutation | 3-5 hours/PC level | 25 GB | 8 GB |

---

## Methods

### Gene Expression Quantification

**Alignment:**
- Tool: minimap2
- Mode: splice-aware (`-ax splice`)
- Strand: forward only (`-uf`)

**Quantification:**
- Tool: featureCounts
- Features: protein-coding and lincRNA genes
- Method: fractional counting for multi-mappers
- Normalization: RPKM

**Filtering:**
- Threshold: expressed in ≥50% of samples
- Result: ~15,000 genes retained

### Transcript Isoform Quantification

**Discovery:**
- Primary method: FLAIR (Full-Length Alternative Isoform analysis of RNA)
- Steps: correct → collapse → quantify
- Minimum support: 10 reads per isoform

**Validation:**
- Tool: SQANTI3
- Purpose: structural classification of novel isoforms
- Categories: FSM, ISM, NIC, NNC

**Filtering:**
- Gene types: protein-coding and lincRNA only
- Expression: ≥50% of samples
- Result: ~40,000 transcripts

### QTL Mapping

**Tool:** QTLtools v1.3

**Approach:**
- Mode: permutation (1000 permutations)
- Window: ±1 Mb from TSS (cis-acting)
- Covariates: expression PCs (optimized)
- FDR control: Benjamini-Hochberg at 5%

**Optimization:**
- Test PC levels: 0, 3, 6, 9, ..., 60
- Select based on: discovery plateau + π₁ statistic
- Typical optimal: 3-15 PCs for gene eQTL, 0-6 for transcript tQTL

**Output:**
- Permutation pass: best SNP per gene/transcript
- Nominal pass: all tested SNP-gene pairs

---

## Documentation

Detailed documentation is available for each analysis stage:

- **[Gene Expression](scripts/01_gene_expression_quantification/README.md)** - Alignment, quantification, filtering
- **[Transcript Quantification](scripts/02_transcript_quantification/README.md)** - Isoform discovery and validation
- **[QTL Mapping](scripts/03_qtl_mapping/README.md)** - eQTL and tQTL discovery workflow
- **[m⁶A Analysis](scripts/04_m6a_modification_analysis/README.md)** - Modification detection and QTL mapping
- **[Replication](scripts/05_replication_analysis/README.md)** - Validation in independent cohorts

---

## Key Features

**Multi-level Analysis**
- Gene expression
- Transcript isoforms
- RNA modifications
- Integrated QTL discovery

**Quality Control**
- Systematic covariate optimization
- Multiple validation methods
- Replication in independent cohort

**Reproducibility**
- Well-documented code
- Complete environment specifications
- Modular pipeline design

**Flexibility**
- Easy parameter adjustment
- Multiple quantification methods
- Adaptable to other datasets

---

## Software Requirements

### Core Dependencies

| Category | Tools |
|----------|-------|
| **Alignment** | minimap2, samtools |
| **Quantification** | featureCounts, FLAIR, Bambu, StringTie |
| **Validation** | SQANTI3 |
| **QTL Mapping** | QTLtools, bgzip, tabix |
| **Statistical Analysis** | R (edgeR, qvalue, data.table) |

### R Packages

```r
# Bioconductor
BiocManager::install(c("edgeR", "rtracklayer", "GenomicRanges", "bambu"))

# CRAN
install.packages(c("data.table", "ggplot2", "dplyr", "tidyr", 
                   "qvalue", "gplots", "corrplot"))
```

### Python Packages

```bash
pip install nanopolish m6anet pandas numpy scipy matplotlib
```

---

## Citation

If you use this pipeline, please cite:

```bibtex
@article{yourpaper2024,
  title={Long-read RNA sequencing reveals transcript-level QTLs in lymphoblastoid cell lines},
  author={Your Name and Collaborators},
  journal={Journal Name},
  year={2024},
  doi={10.xxxx/xxxxx}
}
```

### Software Citations

**Core Tools:**
- minimap2: Li, H. (2018). *Bioinformatics*, 34(18), 3094-3100.
- FLAIR: Tang, A.D. et al. (2020). *Nature Communications*, 11(1), 1438.
- QTLtools: Delaneau, O. et al. (2017). *Nature Communications*, 8, 15452.
- SQANTI3: Tardaguila, M. et al. (2018). *Genome Research*, 28(7), 1096-1108.

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Submit a pull request

For bug reports or questions, please [open an issue](https://github.com/yourusername/DRS-LCLs-qtl-analysis/issues).

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

**Principal Investigator:** Your Name  
**Institution:** Your University  
**Email:** your.email@institution.edu  
**Lab Website:** https://yourlab.org

**For questions:** Use [GitHub Issues](https://github.com/yourusername/DRS-LCLs-qtl-analysis/issues)

---

## Acknowledgments

**Funding:**
- Grant Agency (Award #XXXXX)
- Additional funding sources

**Resources:**
- HPC cluster computational resources
- Sequencing core facilities

**Collaborators:**
- List key contributors

---

<div align="center">

**If you find this work useful, please cite our paper and star this repository**

[Report Bug](https://github.com/yourusername/DRS-LCLs-qtl-analysis/issues) · 
[Request Feature](https://github.com/yourusername/DRS-LCLs-qtl-analysis/issues) · 
[Documentation](https://github.com/yourusername/DRS-LCLs-qtl-analysis/wiki)

</div>
