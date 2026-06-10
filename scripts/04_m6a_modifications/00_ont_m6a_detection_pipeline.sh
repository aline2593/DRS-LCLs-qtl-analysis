#!/usr/bin/env bash
# =============================================================================
# ONT RNA m6A Detection Pipeline
# =============================================================================
# Description: End-to-end pipeline for m6A detection from Oxford Nanopore
#              direct RNA sequencing data.
#
# Steps:
#   1. Transcriptome alignment        (minimap2)
#   2. SAM → BAM conversion          (samtools)
#   3. BAM sorting                   (samtools)
#   4. BAM indexing                  (samtools)
#   4.5 fast5 → slow5 conversion     (slow5tools) [optional, see note below]
#   5. Event alignment               (f5c eventalign)
#   6. m6anet data preparation       (m6anet dataprep)
#   7. m6A inference                 (m6anet inference)
#   8. Merge per-sample outputs      (Python/pandas)
#
# Requirements:
#   - minimap2 v2.17+
#   - samtools v1.23
#   - f5c v1.2
#   - m6anet v2.1.0
#   - slow5tools (required if using slow5/blow5 input for f5c eventalign)
#     Install: https://github.com/hasindu2008/slow5tools
#   - conda environment: ont_env
#
# Note on fast5/slow5:
#   f5c eventalign requires raw signal data in one of two forms:
#     (a) fast5 files indexed with `f5c index` (Option A), OR
#     (b) slow5/blow5 files passed via --slow5 (Option B, recommended).
#   This pipeline uses Option B by default. Run Step 4.5 to convert
#   fast5 → slow5 before running Step 5, unless your basecaller already
#   produced a .index file alongside your fastq (check for ${fastq}.index).
#
# Usage:
#   bash ont_m6a_pipeline.sh            # Submit Steps 1-7
#   bash ont_m6a_pipeline.sh --merge-only  # Run Step 8 after all jobs finish
#
# Configuration: Edit the variables in the CONFIG section below.
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIG
# =============================================================================

FASTQ_DIR="/path/to/fastq_file"
FASTQ_PATTERN="${FASTQ_DIR}/*_all_noduplicate.fastq"

REFERENCE_CDNA="Homo_sapiens.GRCh37.67.cdna.all.fa"

MINIMAP2="/path/to/minimap2-2.17_x64-linux/minimap2"
SAMTOOLS="${HOME}/bin/samtools"

# fast5 / slow5 paths (for f5c eventalign — see Step 4.5 and Step 5)
FAST5_DIR="/path/to/fast5_dir"   # original fast5 directory
SLOW5_DIR="/path/to/slow5_dir"   # output directory for converted slow5/blow5 files

EVENTALIGN_OUT="/path/to/output/eventalign"
M6ANET_OUT="/path/to/output/m6anet"
MERGED_OUT="${M6ANET_OUT}/Nanopore_ModificationsRatio.txt"

SLURM_PARTITION="public-cpu"
SLURM_MEM="10G"
# NOTE: --time limits are not set. Add --time=HH:MM:SS to each sbatch call
#       according to your cluster policy.
#       Typical ranges: 2h (Steps 1-4), 4-8h (Step 4.5), 24-48h (Steps 5-7).

# =============================================================================
# HELPERS
# =============================================================================

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

get_sample_id() {
    # Extracts sample ID from filename (everything before the first underscore)
    local filepath="$1"
    basename "$filepath" | sed 's/_.*//'
}

# =============================================================================
# STEP 1: Transcriptome Alignment
# =============================================================================

log "Step 1: Transcriptome alignment with minimap2"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting alignment for: ${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="align_${sample}" \
        --wrap="
            echo 'Aligning sample: ${sample}'
            ${MINIMAP2} \
                -ax map-ont \
                -t 10 \
                ${REFERENCE_CDNA} \
                ${fastq} \
                > alignment_transcriptome_${sample}.sam
        "
done

# =============================================================================
# STEP 2: SAM → BAM Conversion
# =============================================================================

log "Step 2: SAM to BAM conversion"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting SAM→BAM for: ${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="sam2bam_${sample}" \
        --wrap="
            echo 'Converting SAM to BAM: ${sample}'
            samtools view -S -b \
                alignment_transcriptome_${sample}.sam \
                > alignment_transcriptome_${sample}.bam
        "
done

# =============================================================================
# STEP 3: BAM Sorting
# =============================================================================

log "Step 3: BAM sorting"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting sort for: ${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="sort_${sample}" \
        --wrap="
            echo 'Sorting BAM: ${sample}'
            ${SAMTOOLS} sort \
                -o sorted.alignment_transcriptome_${sample}.bam \
                alignment_transcriptome_${sample}.bam
        "
done

# =============================================================================
# STEP 4: BAM Indexing
# =============================================================================

log "Step 4: BAM indexing"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting index for: ${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="index_${sample}" \
        --wrap="
            echo 'Indexing BAM: ${sample}'
            samtools index sorted.alignment_transcriptome_${sample}.bam
        "
done

# =============================================================================
# STEP 4.5: Raw signal preparation for f5c eventalign
# =============================================================================
# f5c eventalign requires access to raw signal data. Choose ONE option:
#
#   Option A — f5c index (fast5 files):
#     Index your fast5 files directly with f5c index.
#     Use this if your basecaller produced fast5 files and you do NOT
#     want to convert to slow5. Check first whether ${fastq}.index already
#     exists — some basecallers (e.g. Guppy) produce it automatically.
#     If it exists, skip this step entirely.
#     In Step 5: remove --slow5 and add -d ${FAST5_DIR}
#
#   Option B — slow5/blow5 conversion (recommended):
#     Convert fast5 → slow5/blow5 using slow5tools. More portable and
#     faster for large datasets. Only needs to be run once per dataset.
#     Install slow5tools: https://github.com/hasindu2008/slow5tools
#     In Step 5: keep --slow5 ${SLOW5_DIR}/${sample}.blow5
#
# Uncomment the option you want to use. Only run ONE of them.
# Wait for this job to complete before submitting Step 5.
# =============================================================================

# ── Option A: f5c index (fast5) ──────────────────────────────────────────────
# log "Step 4.5 (Option A): Indexing fast5 files with f5c index"
#
# for fastq in ${FASTQ_PATTERN}; do
#     sample=$(get_sample_id "$fastq")
#     log "  Submitting f5c index for: ${sample}"
#
#     sbatch \
#         --partition="${SLURM_PARTITION}" \
#         --mem="${SLURM_MEM}" \
#         --job-name="f5c_index_${sample}" \
#         --wrap="
#             echo 'Indexing fast5 for sample: ${sample}'
#             conda activate ont_env
#             f5c index \
#                 -d ${FAST5_DIR} \
#                 ${fastq}
#             echo 'Indexing complete for: ${sample}'
#         "
# done

# ── Option B: slow5/blow5 conversion (default) ───────────────────────────────
log "Step 4.5 (Option B): Converting fast5 to slow5/blow5 (slow5tools)"
mkdir -p "${SLOW5_DIR}"

sbatch \
    --partition="${SLURM_PARTITION}" \
    --mem="${SLURM_MEM}" \
    --job-name="fast5_to_slow5" \
    --wrap="
        echo 'Converting fast5 → slow5/blow5'
        slow5tools f2s ${FAST5_DIR} \
            -d ${SLOW5_DIR} \
            --iop 8
        echo 'Conversion complete. Output: ${SLOW5_DIR}'
    "

# =============================================================================
# STEP 5: Event Alignment (f5c eventalign)
# =============================================================================
# By default uses slow5/blow5 via --slow5 (Option B from Step 4.5).
#
# If you used Option A (f5c index), replace:
#   --slow5 ${SLOW5_DIR}/${sample}.blow5
# with:
#   -d ${FAST5_DIR}
#
# NOTE: Sorted BAM files (Steps 3-4) AND raw signal data (Step 4.5)
#       must both be ready before submitting this step.
# =============================================================================

log "Step 5: f5c eventalign"
mkdir -p "${EVENTALIGN_OUT}"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting eventalign for: ${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="eventalign_${sample}" \
        --wrap="
            echo 'Running eventalign: ${sample}'
            conda activate ont_env
            f5c eventalign \
                -b sorted.alignment_transcriptome_${sample}.bam \
                -g ${REFERENCE_CDNA} \
                -r ${fastq} \
                --slow5 ${SLOW5_DIR}/${sample}.blow5 \
                --scale-events \
                --iop 50 \
                --rna \
                --signal-index \
                --summary ${EVENTALIGN_OUT}/${sample}_summary.txt \
                --threads 50 \
                > ${EVENTALIGN_OUT}/${sample}_eventalign.tsv
        "
done

# =============================================================================
# STEP 6: m6anet Data Preparation
# =============================================================================

log "Step 6: m6anet dataprep"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting m6anet dataprep for: ${sample}"

    mkdir -p "${M6ANET_OUT}/${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="m6anet_prep_${sample}" \
        --wrap="
            echo 'Running m6anet dataprep: ${sample}'
            m6anet dataprep \
                --eventalign ${EVENTALIGN_OUT}/${sample}_eventalign.tsv \
                --out_dir ${M6ANET_OUT}/${sample} \
                --n_processes 4
        "
done

# =============================================================================
# STEP 7: m6A Inference
# =============================================================================

log "Step 7: m6anet inference"

for fastq in ${FASTQ_PATTERN}; do
    sample=$(get_sample_id "$fastq")
    log "  Submitting m6anet inference for: ${sample}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --mem="${SLURM_MEM}" \
        --job-name="m6anet_inf_${sample}" \
        --wrap="
            echo 'Running m6anet inference: ${sample}'
            m6anet inference \
                --input_dir ${M6ANET_OUT}/${sample} \
                --out_dir ${M6ANET_OUT}/${sample} \
                --n_processes 4 \
                --num_iterations 1000
        "
done

log "All jobs submitted successfully."
log ""
log "REMINDER: Once all SLURM inference jobs (Step 7) have finished, run:"
log "   bash $(basename "$0") --merge-only"
log "This will merge per-sample outputs into: ${MERGED_OUT}"
log "Then run the R filtering script: Rscript m6a_filter_and_coords.R"

# =============================================================================
# STEP 8: Merge Per-Sample m6anet Outputs into a Single Matrix
# =============================================================================
# NOTE: This step must only be run AFTER all Step 7 inference jobs have
#       completed. It collects one data.site_proba.csv per sample and merges
#       them into a single wide matrix (one row per site, one column per sample)
#       ready for the downstream R filtering script.
#
# Run manually once all SLURM jobs are done:
#   bash ont_m6a_pipeline.sh --merge-only
# =============================================================================

merge_m6anet_outputs() {
    log "Step 8: Merging per-sample m6anet outputs..."

    # Verify at least one inference output exists before proceeding
    if ! ls "${M6ANET_OUT}"/*/data.site_proba.csv 1>/dev/null 2>&1; then
        log "  ERROR: No data.site_proba.csv files found under ${M6ANET_OUT}/"
        log "         Make sure all Step 7 inference jobs have completed."
        exit 1
    fi

    python3 - <<EOF
import pandas as pd
import glob
import os

m6anet_out = "${M6ANET_OUT}"
merged_out = "${MERGED_OUT}"

files = sorted(glob.glob(os.path.join(m6anet_out, "*", "data.site_proba.csv")))

if not files:
    raise FileNotFoundError(f"No data.site_proba.csv files found under {m6anet_out}/")

print(f"  Found {len(files)} sample(s) to merge:")
frames = []
for f in files:
    sample = os.path.basename(os.path.dirname(f))
    print(f"    - {sample}")
    df = pd.read_csv(f)

    # Rename the probability column to the sample ID so each sample
    # becomes its own column in the merged matrix
    df = df.rename(columns={"probability_modified": sample})

    # Index on the site identifier columns shared across all samples
    index_cols = ["transcript_id", "transcript_position"]
    frames.append(df.set_index(index_cols))

merged = pd.concat(frames, axis=1).reset_index()

# Reorder: index columns first, then one probability column per sample
merged.to_csv(merged_out, sep="\t", index=False)
print(f"  Merged {len(files)} samples -> {len(merged)} sites total")
print(f"  Output written: {merged_out}")
EOF

    log "Step 8 complete. Merged file: ${MERGED_OUT}"
    log "Next step: run m6a_filter_and_coords.R using ${MERGED_OUT} as input."
}

# Run merge only if --merge-only flag is passed
if [[ "${1:-}" == "--merge-only" ]]; then
    merge_m6anet_outputs
fi
