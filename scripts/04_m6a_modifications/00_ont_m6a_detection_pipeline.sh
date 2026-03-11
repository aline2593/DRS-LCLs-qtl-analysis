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
#   5. Event alignment               (f5c eventalign)
#   6. m6anet data preparation       (m6anet dataprep)
#   7. m6A inference                 (m6anet inference)
#
# Requirements:
#   - minimap2 v2.17+
#   - samtools v1.23
#   - f5c v1.2
#   - m6anet v2.1.0
#   - conda environment: ont_env
#
# Usage:
#   bash ont_m6a_pipeline.sh
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

EVENTALIGN_OUT="/path/to/output/eventalign"
M6ANET_OUT="/path/to/output/m6anet"

SLURM_PARTITION="public-cpu"
SLURM_MEM="10G"          # Increase for Steps 5-7 (eventalign/m6anet); cluster-dependent
# NOTE: --time limits are not set. Add --time=HH:MM:SS to each sbatch call
#       according to your cluster policy. Typical ranges: 2h (Steps 1-4), 24-48h (Steps 5-7).

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
# STEP 5: Event Alignment (f5c eventalign)
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
                -b /path/to/sorted/bam/file/sorted.alignment_transcriptome_${sample}.bam \
                -g ${REFERENCE_CDNA} \
                -r ${fastq} \
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

# =============================================================================
# STEP 8: Merge Per-Sample m6anet Outputs into a Single Matrix
# =============================================================================
# NOTE: This step must only be run AFTER all Step 7 inference jobs have
#       completed. It collects one data.site_ratio.csv per sample and merges
#       them into a single wide matrix (one row per site, one column per sample)
#       ready for the downstream R filtering script.
#
# Run manually once all SLURM jobs are done:
#   bash ont_m6a_pipeline.sh --merge-only
# Or call merge_m6anet_outputs() directly in an interactive session.
# =============================================================================

MERGED_OUT="${M6ANET_OUT}/Nanopore_ModificationsRatio.txt"

merge_m6anet_outputs() {
    log "Step 8: Merging per-sample m6anet outputs..."

    # Verify at least one inference output exists before proceeding
    if ! ls "${M6ANET_OUT}"/*/data.site_ratio.csv 1>/dev/null 2>&1; then
        log "  ERROR: No data.site_ratio.csv files found under ${M6ANET_OUT}/"
        log "         Make sure all Step 7 inference jobs have completed."
        exit 1
    fi

    python3 - <<EOF
import pandas as pd
import glob
import os

m6anet_out = "${M6ANET_OUT}"
merged_out = "${MERGED_OUT}"

files = sorted(glob.glob(os.path.join(m6anet_out, "*", "data.site_ratio.csv")))

if not files:
    raise FileNotFoundError(f"No data.site_ratio.csv files found under {m6anet_out}/")

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

# Run merge automatically only if --merge-only flag is passed,
# otherwise print a reminder so the user knows to run it after jobs finish.
if [[ "${1:-}" == "--merge-only" ]]; then
    merge_m6anet_outputs
else
    log ""
    log "REMINDER: Once all SLURM inference jobs (Step 7) have finished, run:"
    log "   bash $(basename "$0") --merge-only"
    log "This will merge per-sample outputs into: ${MERGED_OUT}"
    log "Then run the R filtering script: Rscript m6a_filter_and_coords.R"
fi
