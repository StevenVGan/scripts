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

r1_files=( "${RAW_DIR}"/*_R1.fastq.gz )
if (( ${#r1_files[@]} == 0 )); then
  echo "[ERROR] No *_R1.fastq.gz files found in $RAW_DIR"
  exit 1
fi

for R1 in "${r1_files[@]}"; do
  R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
  if [[ ! -f "$R2" ]]; then
    echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
    continue
  fi

  sample_name="$(basename "$R1" _R1.fastq.gz)"
  echo "[STEP1] Trimming: $sample_name"

  trim_galore \
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

echo "[STEP1] FastQC on trimmed reads..."
mkdir -p "${TRIM_DIR}/fastqc"

trimmed=( "${TRIM_DIR}"/*_val_*.fq.gz )
if (( ${#trimmed[@]} == 0 )); then
  echo "[ERROR] No trimmed FASTQs (*_val_*.fq.gz) found in $TRIM_DIR"
  exit 1
fi

fastqc \
  --outdir "${TRIM_DIR}/fastqc" \
  --threads "$TRIM_CPU" \
  "${trimmed[@]}"

echo "=== STEP 1 complete ==="
