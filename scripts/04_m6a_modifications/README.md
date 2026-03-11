# m⁶A RNA Modification QTL Pipeline

Detection and genetic mapping of N6-methyladenosine (m⁶A) RNA modifications from Oxford Nanopore direct RNA sequencing (DRS) data, followed by replication in independent short-read and m⁶A-seq datasets.

---

## Pipeline overview

```
00_ont_m6a_detection_pipeline.sh
        ↓  per-sample SLURM/LSF jobs → merge
01_m6a_filter_and_coords.R
        ↓  filtered BED with genomic coordinates
02_maQTL_mapping_pipeline.R
        ↓  BED prep → PCA → permutation pass → FDR → nominal pass
03_maQTL_replication_QTL.R
        ↓  eQTL + sQTL replication (short-read datasets)
04_maQTL_yoruba_replication.R
        ↓  Yoruba m6A-seq peak overlap + m6A-QTL replication
```

---

## Scripts

### `00_ont_m6a_detection_pipeline.sh`
End-to-end m⁶A detection from raw FASTQ files.  
Runs minimap2 alignment → SAM/BAM conversion → f5c eventalign → m6anet dataprep → m6anet inference, each step submitted as a SLURM job per sample.  
After all inference jobs complete, run with `--merge-only` to produce the merged modification ratio matrix.

```bash
# Submit all jobs
bash 00_ont_m6a_detection_pipeline.sh

# After all jobs finish — merge outputs
bash 00_ont_m6a_detection_pipeline.sh --merge-only
```

**Output:** `Nanopore_ModificationsRatio.txt`

---

### `01_m6a_filter_and_coords.R`
Filters modification sites and maps transcript-relative positions to genomic coordinates.  
Retains modifications with ≥50% sample coverage and probability >0.9, restricts to protein-coding genes and lncRNAs (GENCODE v46), and performs exon-aware transcript-to-genome coordinate conversion using cumulative exon lengths.

```bash
Rscript 01_m6a_filter_and_coords.R
```

**Input:** `Nanopore_ModificationsRatio.txt`, `filtered_gencode.v46.annotation.txt`  
**Output:** `filtered_m6a_modifications.txt`, `Motifs_position.bed`, `Ref_file_readytouse.txt`

---

### `02_maQTL_mapping_pipeline.R`
Full m⁶A-QTL mapping pipeline run in sequential stages via flags.

```bash
Rscript 02_maQTL_mapping_pipeline.R                 # BED prep + PCA + submit permutation jobs
Rscript 02_maQTL_mapping_pipeline.R --merge-perm    # merge permutation chunks
Rscript 02_maQTL_mapping_pipeline.R --fdr           # FDR calibration + PC selection plots
Rscript 02_maQTL_mapping_pipeline.R --nominal       # submit nominal pass jobs
Rscript 02_maQTL_mapping_pipeline.R --merge-nominal # merge nominal chunks + annotate gene IDs
```

**Input:** `Nanopore_ModificationsRatio.txt`, `merged.pcs.0`
**Output:** `Nanopore_ModificationsRatio_NEW.bed.gz`, `genes.50percent.pca`, `merged.pcs`, `m6A_perm1000_cov{N}_merged.txt.gz`, `pval_dists.pdf`, `no_sig_by_pcs_FDRs.pdf`, `LCLs_m6A_FDR{5,10}_{N}PCs.significant.txt`, `Nanopore_ModificationsRatio_nominal.bed.gz`, `Nominal_gene_ID.txt`

---

### `03_maQTL_replication_QTL.R`
Replication of significant m⁶A-QTLs in two independent short-read QTL datasets.  
Matches maQTL nominal p-values at overlapping gene–SNP pairs, applies q-value FDR correction, and reports significant replicating hits.

- **Illumina eQTL** — 7,658 significant eQTLs from Delaneau et al., Science 2019 => `LCL_RNA.chunkALL.significant_permutation.txt.gz`
- **GTEx sQTL** — splicing QTLs from EBV-transformed lymphocytes => `Cells_EBV-transformed_lymphocytes.v8.sgenes_mod.txt`

```bash
Rscript 03_maQTL_replication_QTL.R
```

**Input:** `Nominal_gene_ID.txt`, `LCLs_m6A_FDR10_3PCs.significant.txt`, `LCL_RNA.chunkALL.significant_permutation.txt.gz`, `Cells_EBV-transformed_lymphocytes.v8.sgenes_mod.txt`
**Output:** `FDR5_eQTLs_gene_recap_modifications.txt`, `FDR5_sQTLs_gene_recap_modifications.txt`, `QTL_replication_summary.txt`, 

---

### `04_maQTL_yoruba_replication.R`
Orthogonal validation using Yoruba LCL m⁶A-seq data (1000 Genomes Project).  
Overlaps DRS modification sites with Yoruba m⁶A-seq peaks (by position and by position + gene), computes overlap rates across three modification sets (all / m6aQTL-tested / significant m6aQTL hits), checks replication of lead SNP–gene pairs in the Yoruba m⁶A-QTL dataset and describe the CALR gene example.

```bash
Rscript 04_maQTL_yoruba_replication.R
```

**Input:** `m6AQTL.m6APeak_logOR_GC.IP.adjusted_qqnorm.15PCs.fastQTL.nominals.rds`, `m6APeak_logOR_GC.IP.adjusted_qqnorm.fastQTL.txt.gz`, `filtered_m6a_modifications.txt`, `LCLs_m6A_FDR10_3PCs.significant.txt`
**Output:** `overlap_summary_stats.txt`, `nano_overlap_peaks_by_position.txt`, `nano_overlap_peaks_by_position_AND_gene.txt`,`maQTL_hits_yoruba_peak_status.txt`, `CALR_replication_yoruba.txt`, `upset_modifications_overlap_maQTL_tested_sig.pdf`

---

## Dependencies

| Tool | Version | Used in |
|------|---------|---------|
| minimap2 | ≥2.17 | 00 |
| samtools | any | 00 |
| f5c | any | 00 |
| m6anet | any | 00 |
| QTLtools | 1.3 | 02 |
| bcftools / tabix / bgzip | any | 02 |
| R ≥4.1 | — | 01–04 |
| data.table, ggplot2, qvalue | — | 01–04 |
| GenomicRanges, IRanges | — | 01, 04 |
| UpSetR | — | 04 |

---
