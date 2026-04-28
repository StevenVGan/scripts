#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 2_bowtie2.sh — PRO-seq alignment + strand-aware bigWig generation
#
# Forked from csRNA/2_bowtie2.sh. The only change is in run_bamcoverage_tracks():
# when PROSEQ_FLIP_STRAND=1, the --filterRNAstrand labels are SWAPPED so that
# ${sample}_fwd.bw represents the nascent RNA on the + strand (PRO-seq R1 is
# the reverse complement of the nascent RNA; see 0_config.sh notes).
#
# Inputs
#   - Trimmed FASTQs in ${TRIM_DIR} (PE: *_R1_val_1.fq.gz + *_R2_val_2.fq.gz; SE: *_trimmed.fq.gz)
#   - Bowtie2 genome index: ${GENOME_INDEX}
#
# Outputs
#   - ${BAM_DIR}/*_sorted.bam + .bai, *_bowtie2.log
#   - ${BAMQC_DIR}/*.idxstats.txt, *.stats
#   - ${TRACK_DIR}/${sample}.bw (optional unstranded, when COMBINED_BIGWIG=1)
#   - ${TRACK_DIR}/${sample}_fwd.bw and _rev.bw (when STRAND_BIGWIG=1; RNA-strand when PROSEQ_FLIP_STRAND=1)
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "2_bowtie2"

check_cmd bowtie2
check_cmd samtools
check_cmd bamCoverage

echo "=== STEP 2: Bowtie2 -> sorted BAM + bigWig + basic BAM QC ==="
echo "[INFO] GENOME_INDEX:       $GENOME_INDEX"
echo "[INFO] TRIM_DIR:           $TRIM_DIR"
echo "[INFO] BAM_DIR:            $BAM_DIR"
echo "[INFO] TRACK_DIR:          $TRACK_DIR"
echo "[INFO] Mode:               $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"
echo "[INFO] BT2_EXTRA:          ${BT2_EXTRA:-}"
echo "[INFO] STRAND_BIGWIG:      ${STRAND_BIGWIG:-0}  COMBINED_BIGWIG: ${COMBINED_BIGWIG:-1}"
echo "[INFO] PROSEQ_FLIP_STRAND: ${PROSEQ_FLIP_STRAND:-0} (swaps fwd/rev bigWig labels so they reflect RNA strand)"

# bigWig(s): optional unstranded + strand-separated.
# When PROSEQ_FLIP_STRAND=1, we swap --filterRNAstrand labels so the output
# file named "_fwd.bw" contains reads from the RNA's + strand (which, for
# PRO-seq, are the reads that aligned to the - strand).
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

  local strand_for_fwd strand_for_rev
  if [[ "${PROSEQ_FLIP_STRAND:-0}" -eq 1 ]]; then
    # PRO-seq: R1 is antisense of nascent RNA → swap filterRNAstrand labels
    strand_for_fwd="reverse"
    strand_for_rev="forward"
  else
    strand_for_fwd="forward"
    strand_for_rev="reverse"
  fi

  if [[ ! -f "$bw_fwd" ]]; then
    bamCoverage "${common[@]}" --filterRNAstrand "$strand_for_fwd" -o "$bw_fwd"
  fi
  if [[ ! -f "$bw_rev" ]]; then
    bamCoverage "${common[@]}" --filterRNAstrand "$strand_for_rev" -o "$bw_rev"
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
