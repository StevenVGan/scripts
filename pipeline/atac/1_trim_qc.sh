#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 1_trim_qc.sh — ATAC adapter/quality trimming + FastQC
#
# Difference vs cutrun: passes $TRIM_ADAPTER_PRESET (defaults to --nextera)
# so Trim Galore looks for the Tn5 ME adapter (CTGTCTCTTATA) instead of TruSeq.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "1_trim_qc"

check_cmd trim_galore
check_cmd cutadapt
check_cmd fastqc

echo "=== STEP 1: Trimming + FastQC ==="
echo "[INFO] RAW_DIR:  $RAW_DIR"
echo "[INFO] TRIM_DIR: $TRIM_DIR"
echo "[INFO] Adapter:  $TRIM_ADAPTER_PRESET"
echo "[INFO] Mode:     $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"

r1_files=( "${RAW_DIR}"/*_R1.fastq.gz )
if (( ${#r1_files[@]} == 0 )); then
  echo "[ERROR] No *_R1.fastq.gz files found in $RAW_DIR"
  exit 1
fi

# shellcheck disable=SC2206
adapter_args=( $TRIM_ADAPTER_PRESET )

if [[ "${SE:-0}" -eq 1 ]]; then
  for R1 in "${r1_files[@]}"; do
    sample_name="$(basename "$R1" _R1.fastq.gz)"
    echo "[STEP1] Trimming (single-end): $sample_name"
    trim_galore \
      "${adapter_args[@]}" \
      --cores "$TRIM_CPU" \
      --output_dir "$TRIM_DIR" \
      -q "$TRIM_QUAL" \
      --phred33 \
      --stringency "$TRIM_STRINGENCY" \
      --length "$TRIM_MIN_LENGTH" \
      --gzip \
      "$R1"
  done
  trimmed=( "${TRIM_DIR}"/*_trimmed.fq.gz )
else
  for R1 in "${r1_files[@]}"; do
    R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
      continue
    fi
    sample_name="$(basename "$R1" _R1.fastq.gz)"
    echo "[STEP1] Trimming: $sample_name"
    trim_galore \
      "${adapter_args[@]}" \
      --paired \
      --cores "$TRIM_CPU" \
      --output_dir "$TRIM_DIR" \
      -q "$TRIM_QUAL" \
      --phred33 \
      --stringency "$TRIM_STRINGENCY" \
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
  --threads "$TRIM_CPU" \
  "${trimmed[@]}"

echo "=== STEP 1 complete ==="
