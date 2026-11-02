#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 1_trim_qc.sh
#
# Purpose
#   Step 1 of the pipeline: adapter/quality trimming and FastQC on trimmed reads.
#   MultiQC is run in step 5 (5_qc.sh).
#
# Inputs
#   - Raw paired FASTQs in:
#       ${RAW_DIR}/*_R1.fastq.gz and ${RAW_DIR}/*_R2.fastq.gz
#
# Outputs
#   - Trim Galore paired outputs in ${TRIM_DIR}:
#       *_R1_val_1.fq.gz and *_R2_val_2.fq.gz
#   - FastQC reports in:
#       ${TRIM_DIR}/fastqc/
#
# Behavior
#   - Loops over *_R1.fastq.gz, requires the paired *_R2.fastq.gz.
#   - Runs Trim Galore with parameters from 0_config.sh.
#   - Stringency per library: TRIM_STRINGENCY for paired-end libraries whose shorter mate is at least
#     TRIM_STRINGENCY_MIN_READ bp (longest of the first 1,000 reads); TRIM_STRINGENCY_FALLBACK for
#     single-end runs and shorter reads. The choice is logged per library.
#   - Runs FastQC on all *_val_*.fq.gz outputs.
#
# Logging
#   - Writes a timestamped log to ${LOG_DIR}/1_trim_qc_*.log
#
# Requirements
#   - trim_galore, cutadapt, fastqc on PATH.
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "1_trim_qc"

check_cmd trim_galore
check_cmd cutadapt
check_cmd fastqc

echo "=== STEP 1: Trimming + FastQC (no MultiQC) ==="
echo "[INFO] RAW_DIR:  $RAW_DIR"
echo "[INFO] TRIM_DIR: $TRIM_DIR"
echo "[INFO] Mode:     $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"

# quality trimming: 2-colour aware (--2colour) or plain -q — see TRIM_2COLOUR in 0_config.sh
if [[ "${TRIM_2COLOUR:-1}" -eq 1 ]]; then TRIM_QARGS=(--2colour "$TRIM_QUAL"); else TRIM_QARGS=(-q "$TRIM_QUAL"); fi
echo "[INFO] Quality:  ${TRIM_QARGS[*]}"
r1_files=( "${RAW_DIR}"/*_R1.fastq.gz )
if (( ${#r1_files[@]} == 0 )); then
  echo "[ERROR] No *_R1.fastq.gz files found in $RAW_DIR"
  exit 1
fi

# Stringency guard (2026-11-01; see TRIM_STRINGENCY in 0_config.sh). Defaults cover a project 0_config.sh from before then.
STRINGENCY_FALLBACK="${TRIM_STRINGENCY_FALLBACK:-5}"
STRINGENCY_MIN_READ="${TRIM_STRINGENCY_MIN_READ:-50}"
# Longest of the first 1,000 reads of a FASTQ(.gz), 0 if it has none. zcat's status is dropped on purpose: awk exits
# early, so zcat dies of SIGPIPE, which pipefail would otherwise turn into a silent exit 141.
fastq_read_len() {
  { zcat -f "$1" 2>/dev/null || true; } | awk 'NR % 4 == 2 && length($0) > m {m = length($0)} NR >= 4000 {exit} END {print m + 0}'
}

if [[ "${SE:-0}" -eq 1 ]]; then
  # Single-end: no R2, output *_trimmed.fq.gz
  for R1 in "${r1_files[@]}"; do
    sample_name="$(basename "$R1" _R1.fastq.gz)"
    stringency="$STRINGENCY_FALLBACK"
    (( TRIM_STRINGENCY > stringency )) && stringency="$TRIM_STRINGENCY"   # never looser than the project's own setting
    echo "[STEP1] Trimming (single-end): $sample_name, stringency $stringency (single-end runs use TRIM_STRINGENCY_FALLBACK)"
    trim_galore \
      --cores "$TRIM_CPU" \
      --output_dir "$TRIM_DIR" \
      "${TRIM_QARGS[@]}" \
      --phred33 \
      --stringency "$stringency" \
      --length "$TRIM_MIN_LENGTH" \
      --gzip \
      "$R1"
  done
  trimmed=( "${TRIM_DIR}"/*_trimmed.fq.gz )
else
  # Paired-end
  for R1 in "${r1_files[@]}"; do
    R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
      continue
    fi
    sample_name="$(basename "$R1" _R1.fastq.gz)"
    len1="$(fastq_read_len "$R1")"; len2="$(fastq_read_len "$R2")"
    read_len=$(( len1 < len2 ? len1 : len2 ))
    stringency="$TRIM_STRINGENCY"
    if (( read_len < STRINGENCY_MIN_READ && STRINGENCY_FALLBACK > stringency )); then
      stringency="$STRINGENCY_FALLBACK"
      echo "[WARN] $sample_name: reads ${len1}/${len2} bp, under ${STRINGENCY_MIN_READ} bp: stringency $stringency instead of $TRIM_STRINGENCY." \
           "Expect a hole in the fragment lengths at read length - 4 ... - 1 (5_qc.sh frag_hole_check)."
    fi
    echo "[STEP1] Trimming: $sample_name, reads ${len1}/${len2} bp, stringency $stringency"
    trim_galore \
      --paired \
      --cores "$TRIM_CPU" \
      --output_dir "$TRIM_DIR" \
      "${TRIM_QARGS[@]}" \
      --phred33 \
      --stringency "$stringency" \
      --length "$TRIM_MIN_LENGTH" \
      --gzip \
      "$R1" "$R2"
  done
  trimmed=( "${TRIM_DIR}"/*_val_*.fq.gz )
fi

echo "[STEP1] FastQC on trimmed reads..."
mkdir -p "${TRIM_DIR}/fastqc"

if (( ${#trimmed[@]} == 0 )); then
  echo "[ERROR] No trimmed FASTQs found in $TRIM_DIR"
  exit 1
fi

fastqc \
  --outdir "${TRIM_DIR}/fastqc" \
  --threads "${FASTQC_THREADS:-$TRIM_CPU}" \
  "${trimmed[@]}"

echo "=== STEP 1 complete ==="
