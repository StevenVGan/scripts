#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 1_trim_qc.sh — Trim Galore + FastQC
#
# Single Trim Galore call per sample:
#   --nextseq ${TRIM_QUAL}          polyG-aware quality trim (NovaSeq/NextSeq
#                                   2-color chemistry emits G for "no signal")
#   (auto-detect on R1)             Illumina Universal 3′ adapter (NEBNext/TruSeq)
#   --adapter2 ${SMALL_RNA_5P_ADAPTER}
#                                   R2 3′ = Illumina Small RNA 5′ adapter
#                                   readthrough (inserts < read length)
#
# FastQC runs on the final trimmed FASTQs; MultiQC is deferred to 5_qc.sh so
# there is only one aggregate report per project.
#
# Outputs:
#   ${TRIM_DIR}/<sample>_R{1,2}_val_{1,2}.fq.gz             (PE)
#   ${TRIM_DIR}/<sample>_trimmed.fq.gz                      (SE)
#   ${TRIM_DIR}/*_trimming_report.txt
#   ${TRIM_DIR}/fastqc/
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "1_trim_qc"

check_cmd trim_galore
check_cmd cutadapt
check_cmd fastqc

echo "=== STEP 1: Trim Galore + FastQC ==="
echo "[INFO] RAW_DIR:          $RAW_DIR"
echo "[INFO] TRIM_DIR:         $TRIM_DIR"
echo "[INFO] Mode:             $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"
echo "[INFO] --nextseq cutoff: $TRIM_QUAL (NovaSeq/NextSeq polyG-aware quality trim)"
if [[ "${SE:-0}" -ne 1 ]]; then
  echo "[INFO] --adapter2 (R2):  ${SMALL_RNA_5P_ADAPTER:-(disabled, auto-detect both mates)}"
fi

tg_extra=()
if [[ -n "${TRIM_GALORE_EXTRA:-}" ]]; then
  # shellcheck disable=SC2206
  read -r -a tg_extra <<< "$TRIM_GALORE_EXTRA"
fi

r1_files=( "${RAW_DIR}"/*_R1.fastq.gz )
if (( ${#r1_files[@]} == 0 )); then
  echo "[ERROR] No *_R1.fastq.gz files found in $RAW_DIR"
  exit 1
fi

mkdir -p "$TRIM_DIR"

if [[ "${SE:-0}" -eq 1 ]]; then
  for R1 in "${r1_files[@]}"; do
    sample="$(basename "$R1" _R1.fastq.gz)"
    echo "[STEP1] Trim Galore (SE): $sample"
    trim_galore \
      --cores "$TRIM_CPU" \
      --output_dir "$TRIM_DIR" \
      --nextseq "$TRIM_QUAL" \
      --phred33 \
      --stringency "$TRIM_STRINGENCY" \
      --length "$TRIM_MIN_LENGTH" \
      --gzip \
      ${tg_extra[@]+"${tg_extra[@]}"} \
      "$R1"
  done
else
  a2_opts=()
  if [[ -n "${SMALL_RNA_5P_ADAPTER:-}" ]]; then
    a2_opts=( --adapter2 "$SMALL_RNA_5P_ADAPTER" )
  fi
  for R1 in "${r1_files[@]}"; do
    R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
      continue
    fi
    sample="$(basename "$R1" _R1.fastq.gz)"
    echo "[STEP1] Trim Galore (PE): $sample"
    trim_galore \
      --paired \
      --cores "$TRIM_CPU" \
      --output_dir "$TRIM_DIR" \
      --nextseq "$TRIM_QUAL" \
      --phred33 \
      --stringency "$TRIM_STRINGENCY" \
      --length "$TRIM_MIN_LENGTH" \
      --gzip \
      "${a2_opts[@]}" \
      ${tg_extra[@]+"${tg_extra[@]}"} \
      "$R1" "$R2"
  done
fi

if [[ "${SE:-0}" -eq 1 ]]; then
  trimmed=( "${TRIM_DIR}"/*_trimmed.fq.gz )
else
  trimmed=( "${TRIM_DIR}"/*_val_*.fq.gz )
fi
if (( ${#trimmed[@]} == 0 )); then
  echo "[ERROR] No trimmed FASTQs found in $TRIM_DIR"
  exit 1
fi

fastqc_dir="${TRIM_DIR}/fastqc"
mkdir -p "$fastqc_dir"
rm -f "${fastqc_dir}"/*.html "${fastqc_dir}"/*.zip

echo "[STEP1] FastQC on trimmed reads..."
fastqc --outdir "$fastqc_dir" --threads "$TRIM_CPU" "${trimmed[@]}"

echo "=== STEP 1 complete ==="
