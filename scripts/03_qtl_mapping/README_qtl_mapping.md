# cis-QTL Mapping (eQTL and tQTL)

This directory contains scripts for mapping expression and transcript quantitative trait loci (QTLs) using QTLtools (https://qtltools.github.io/qtltools/).

## Quick Overview

| Analysis | What it does | Output |
|----------|-------------|--------|
| **eQTL** | Maps genetic variants affecting gene expression | Genes with expression QTLs |
| **tQTL** | Maps genetic variants affecting transcript isoform usage | Transcripts with isoform QTLs |

---

## Scripts

### eQTL Mapping

**1. `01_eqtl_permutation.sh` - eQTL Discovery**
- Maps all gene-level expression QTLs
- Optimizes number of expression PCs as covariates
- Tests PC levels: 0, 3, 6, ..., 60
- Outputs significant eGenes at FDR < 5%

**2. `02_eqtl_nominal.sh` - Complete eQTL Results**
- Gets p-values for ALL gene-SNP pairs
- Required for colocalization and fine-mapping
- Uses optimal PC number from permutation step

### tQTL Mapping

**3. `03_tqtl_permutation.sh` - tQTL Discovery**
- Maps transcript-level isoform QTLs
- Same approach as eQTL but for transcripts
- Finds isoform-switching variants
- Usually needs fewer PCs than eQTL

**4. `04_tqtl_nominal.sh` - Complete tQTL Results**
- Gets p-values for ALL transcript-SNP pairs
- Required for isoform-switching analysis
- Larger output than eQTL (more transcripts)

---

## Workflow

### Complete eQTL Analysis
```bash
# Step 1: PC optimization and discovery
bash 01_eqtl_permutation.sh
# → Check plots in output/pc_optimization/eQTL_PC_optimization.pdf
# → Choose optimal PC number (where eGene discovery plateaus)

# Step 2: Nominal pass with optimal PCs
# Edit script to set OPTIMAL_PCS variable
bash 02_eqtl_nominal.sh
```

### Complete tQTL Analysis
```bash
# Step 1: PC optimization and discovery
bash 03_tqtl_permutation.sh
# → Check plots in output/pc_optimization/tQTL_PC_optimization.pdf
# → Choose optimal PC number (often 0-3 for transcripts)

# Step 2: Nominal pass
# Edit script to set OPTIMAL_PCS variable
bash 04_tqtl_nominal.sh
```

---

## Understanding PC Optimization

### Why optimize PCs?
- **Too few PCs**: Technical variation inflates false positives
- **Too many PCs**: Remove true genetic signal

### How to choose optimal PC number?
1. Look at discovery curve: where does eGene/eTranscript count plateau?
2. Check π₁ statistic: proportion of true positive signals
3. Balance power vs false discovery

### Typical Results
| Dataset | Gene eQTL | Transcript tQTL |
|---------|-----------|-----------------|
| Optimal PCs | 3-10 | 0-5 |
| Reason | More noise in gene counts | Cleaner isoform quantification |

---

## Key Parameters

### Common to both eQTL and tQTL
```bash
--permute 1000    # 1000 permutations for FDR control
--chunk X 20      # Parallel processing (20 chunks)
--normal          # Rank-normal transform phenotypes
--cov FILE        # Covariate file (PCs + technical)
```

### tQTL-specific
```bash
--grp-best        # Test best transcript per gene
                  # (avoids multiple testing across isoforms)
```

---

## Output Files

### Permutation Results
```
Column 1:  Phenotype ID/ Group ID(gene)
Column 2:  Chromosome
Column 5.1:  Transcript ID if --grp-best used
Column 6:  Number of variants tested in cis
Column 8: Variants ID
Column 16: nominal p-value
Column 18: Beta (effect size)
Column 19: Adjusted Empirical p-value
Column 20: Adjusted Beta p-value (from permutations)
```

### Nominal Results
```
Column 1:  Phenotype ID/ Group ID(gene)
Column 2:  Chromosome
Column 5.1:  Transcript ID if --grp-best used
Column 6:  Number of variants tested in cis
Column 8: Variants ID
Column 12: nominal p-value
Column 14: Beta (effect size)
Column 15: Best hit (0 or 1)
```
- ONE line per phenotype-SNP pair tested
- Much larger file (all associations, not just best)
```
---
## Dependencies

- QTLtools (v1.3+)
- htslib (bgzip, tabix)
- R packages: qvalue, data.table
- HPC cluster with job scheduler

---
## Notes

### Cis vs Trans
- **Cis-QTL**: Variant near gene (default ±1Mb)
- **Trans-QTL**: Variant far from gene
- These scripts focus on cis-QTL (more power)

### --grp-best Flag
- Groups transcripts by gene
- Tests only best transcript per gene
- Reduces multiple testing burden

