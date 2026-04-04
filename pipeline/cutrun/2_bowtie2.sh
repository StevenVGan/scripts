#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 2_bowtie2.sh
#
# Purpose
#   Step 2 of the pipeline: alignment and BAM/track generation.
#   Performs Bowtie2 alignment, produces sorted+indexed BAMs, generates bigWig
#   tracks, and writes basic BAM QC stats (idxstats, samtools stats).
#
# Inputs
#   - Trimmed paired FASTQs in ${TRIM_DIR}:
#       *_R1_val_1.fq.gz and *_R2_val_2.fq.gz
#   - Bowtie2 genome index: ${GENOME_INDEX} (configured in 0_config.sh)
#
# Outputs
#   - Sorted BAM + index:
#       ${BAM_DIR}/*_sorted.bam and ${BAM_DIR}/*_sorted.bam.bai
#   - Bowtie2 logs:
#       ${BAM_DIR}/*_bowtie2.log
#   - Per-contig mapping stats:
#       ${BAMQC_DIR}/*.idxstats.txt
#   - samtools stats (for MultiQC):
#       ${BAMQC_DIR}/*.stats
#   - bigWig tracks:
#       ${TRACK_DIR}/*.bw (blacklist-filtered when BLACKLIST exists in 0_config.sh)
#
# Behavior
#   - Streams Bowtie2 output directly into samtools sort (does not write SAM).
#   - Skips alignment if *_sorted.bam + .bai already exist.
#   - Skips bamCoverage if the .bw already exists.
#   - Optional: deletes trimmed FASTQs after alignment (space-saving), controlled by
#       DELETE_TRIMMED_AFTER_ALIGN in 0_config.sh.
#
# Logging
#   - Writes a timestamped log to ${LOG_DIR}/2_bowtie2_*.log
#
# Requirements
#   - bowtie2, samtools, bamCoverage on PATH.
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "2_bowtie2"

check_cmd bowtie2
check_cmd samtools
check_cmd bamCoverage

echo "=== STEP 2: Bowtie2 -> sorted BAM + bigWig + basic BAM QC ==="
echo "[INFO] GENOME_INDEX: $GENOME_INDEX"
echo "[INFO] TRIM_DIR:     $TRIM_DIR"
echo "[INFO] BAM_DIR:      $BAM_DIR"
echo "[INFO] TRACK_DIR:   $TRACK_DIR"
echo "[INFO] Mode:        $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"

if [[ "${SE:-0}" -eq 1 ]]; then
  # Single-end: *_trimmed.fq.gz, bowtie2 -U
  trimmed_files=( "${TRIM_DIR}"/*_trimmed.fq.gz )
  if (( ${#trimmed_files[@]} == 0 )); then
    echo "[ERROR] No *_trimmed.fq.gz files found in $TRIM_DIR"
    exit 1
  fi
  for R1 in "${trimmed_files[@]}"; do
    base="$(basename "$R1" .fq.gz)"
    if [[ "$base" == *_R1_trimmed ]]; then
      sample_name="${base%_R1_trimmed}"
    else
      sample_name="${base%_trimmed}"
    fi
    sorted_bam="${BAM_DIR}/${sample_name}_sorted.bam"
    bw_file="${TRACK_DIR}/${sample_name}.bw"
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
      echo "[STEP2] Found existing BAM+BAI for $sample_name, skipping alignment."
    else
      echo "[STEP2] Aligning (single-end): $sample_name"
      bowtie2 -x "$GENOME_INDEX" -U "$R1" --very-sensitive -p "$BT2_CPU" 2> "$bt2_log" \
        | samtools sort -@ "$BT2_CPU" -o "$sorted_bam" -
      samtools index "$sorted_bam"
    fi

    [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
    [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"

    if [[ ! -f "$bw_file" ]]; then
      bam_cov_args=(-b "$sorted_bam" -o "$bw_file" -p "$BAMCOV_CPU" -bs "$BINSIZE" --effectiveGenomeSize "$GENOMESIZE" --normalizeUsing "$NORMALIZE" --ignoreForNormalization "$IGNORE_CHR" --ignoreDuplicates)
      [[ -f "${BLACKLIST:-}" ]] && bam_cov_args+=( --blackListFileName "$BLACKLIST" )
      bamCoverage "${bam_cov_args[@]}"
    fi
  done
  [[ "$DELETE_TRIMMED_AFTER_ALIGN" -eq 1 ]] && rm -f "${TRIM_DIR}"/*_trimmed.fq.gz || true
else
  # Paired-end
  r1_files=( "${TRIM_DIR}"/*_R1_val_1.fq.gz )
  if (( ${#r1_files[@]} == 0 )); then
    echo "[ERROR] No *_R1_val_1.fq.gz files found in $TRIM_DIR"
    exit 1
  fi

  for R1 in "${r1_files[@]}"; do
    R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
      continue
    fi

    sample_name="$(basename "$R1" _R1_val_1.fq.gz)"
    sorted_bam="${BAM_DIR}/${sample_name}_sorted.bam"
    bw_file="${TRACK_DIR}/${sample_name}.bw"
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
      echo "[STEP2] Found existing BAM+BAI for $sample_name, skipping alignment."
    else
      echo "[STEP2] Aligning: $sample_name"
      bowtie2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" --very-sensitive -p "$BT2_CPU" 2> "$bt2_log" \
        | samtools sort -@ "$BT2_CPU" -o "$sorted_bam" -
      samtools index "$sorted_bam"
    fi

    [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
    [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"

    if [[ ! -f "$bw_file" ]]; then
      bam_cov_args=(-b "$sorted_bam" -o "$bw_file" -p "$BAMCOV_CPU" -bs "$BINSIZE" --effectiveGenomeSize "$GENOMESIZE" --normalizeUsing "$NORMALIZE" --ignoreForNormalization "$IGNORE_CHR" --ignoreDuplicates)
      [[ -f "$BLACKLIST" ]] && bam_cov_args+=( --blackListFileName "$BLACKLIST" )
      bamCoverage "${bam_cov_args[@]}"
    fi
  done

  [[ "$DELETE_TRIMMED_AFTER_ALIGN" -eq 1 ]] && rm -f "${TRIM_DIR}"/*_val_*.fq.gz || true
fi

echo "=== STEP 2 complete ==="
