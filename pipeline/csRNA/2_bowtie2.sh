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
#   - bigWig tracks (see STRAND_BIGWIG / COMBINED_BIGWIG in 0_config.sh):
#       ${TRACK_DIR}/${sample}.bw (optional unstranded)
#       ${TRACK_DIR}/${sample}_fwd.bw and _rev.bw (strand-separated when STRAND_BIGWIG=1)
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
echo "[INFO] BT2_EXTRA:   ${BT2_EXTRA:-}"
echo "[INFO] STRAND_BIGWIG=${STRAND_BIGWIG:-0} COMBINED_BIGWIG=${COMBINED_BIGWIG:-1}"

# bigWig(s): optional unstranded + strand-separated (deeptools --filterRNAstrand)
run_bamcoverage_tracks() {
  local sample_name="$1"
  local sorted_bam="$2"
  local common=( -b "$sorted_bam" -p "$BAMCOV_CPU" -bs "$BINSIZE" --effectiveGenomeSize "$GENOMESIZE" --normalizeUsing "$NORMALIZE" --ignoreForNormalization "$IGNORE_CHR" )
  [[ "${BAMCOV_IGNORE_DUP:-0}" -eq 1 ]] && common+=( --ignoreDuplicates )
  [[ -f "${BLACKLIST:-}" ]] && common+=( --blackListFileName "$BLACKLIST" )

  if [[ "${STRAND_BIGWIG:-0}" -ne 1 ]]; then
    local bw_file="${TRACK_DIR}/${sample_name}.bw"
    if [[ ! -f "$bw_file" ]]; then
      bamCoverage "${common[@]}" -o "$bw_file"
    fi
    return
  fi

  if [[ "${COMBINED_BIGWIG:-1}" -eq 1 ]]; then
    local bw_file="${TRACK_DIR}/${sample_name}.bw"
    if [[ ! -f "$bw_file" ]]; then
      bamCoverage "${common[@]}" -o "$bw_file"
    fi
  fi

  local bw_fwd="${TRACK_DIR}/${sample_name}_fwd.bw"
  local bw_rev="${TRACK_DIR}/${sample_name}_rev.bw"
  if [[ ! -f "$bw_fwd" ]]; then
    bamCoverage "${common[@]}" --filterRNAstrand forward -o "$bw_fwd"
  fi
  if [[ ! -f "$bw_rev" ]]; then
    bamCoverage "${common[@]}" --filterRNAstrand reverse -o "$bw_rev"
  fi
}

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
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
      echo "[STEP2] Found existing BAM+BAI for $sample_name, skipping alignment."
    else
      echo "[STEP2] Aligning (single-end): $sample_name"
      bt2_cmd=( bowtie2 -x "$GENOME_INDEX" -U "$R1" --very-sensitive-local -p "$BT2_CPU" )
      if [[ -n "${BT2_EXTRA:-}" ]]; then
        read -r -a bt2_extra_arr <<< "$BT2_EXTRA"
        bt2_cmd+=( "${bt2_extra_arr[@]}" )
      fi
      "${bt2_cmd[@]}" 2> "$bt2_log" \
        | samtools sort -@ "$BT2_CPU" -o "$sorted_bam" -
      samtools index "$sorted_bam"
    fi

    [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
    [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"

    run_bamcoverage_tracks "$sample_name" "$sorted_bam"
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
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
      echo "[STEP2] Found existing BAM+BAI for $sample_name, skipping alignment."
    else
      echo "[STEP2] Aligning: $sample_name"
      bt2_cmd=( bowtie2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" --very-sensitive-local -p "$BT2_CPU" )
      if [[ -n "${BT2_EXTRA:-}" ]]; then
        read -r -a bt2_extra_arr <<< "$BT2_EXTRA"
        bt2_cmd+=( "${bt2_extra_arr[@]}" )
      fi
      "${bt2_cmd[@]}" 2> "$bt2_log" \
        | samtools sort -@ "$BT2_CPU" -o "$sorted_bam" -
      samtools index "$sorted_bam"
    fi

    [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
    [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"

    run_bamcoverage_tracks "$sample_name" "$sorted_bam"
  done

  [[ "$DELETE_TRIMMED_AFTER_ALIGN" -eq 1 ]] && rm -f "${TRIM_DIR}"/*_val_*.fq.gz || true
fi

echo "=== STEP 2 complete ==="
