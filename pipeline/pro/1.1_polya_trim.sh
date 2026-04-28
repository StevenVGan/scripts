#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 1.1_polya_trim.sh — OPTIONAL poly-A / poly-T removal after step 1
#
# Why this step exists
#   Some PRO-seq libraries (notably early MCF7 runs) carry high poly-A content
#   past adapter trimming — most likely mature mRNA / poly(A)+ carryover
#   through the nuclear run-on. Because PRO-seq R1 is the reverse complement
#   of the nascent RNA:
#     - poly(A) at the 3′ end of the ORIGINAL RNA (mRNA tail / mRNA
#       contamination) → appears as poly-T at the 5′ end of R1;
#     - poly(A) visible in the READS (FastQC "AAAA..." overrepresented) arises
#       when inserts span into genomic A-tracts or when the library has sense-
#       oriented contaminants — trim these from the 3′ end of R1 as well.
#   Trimming both poly-A (3′) and poly-T (5′) is the safe default.
#
# Input
#   Already-trimmed step-1 output:
#     SE: ${TRIM_DIR}/*_trimmed.fq.gz
#     PE: ${TRIM_DIR}/*_R1_val_1.fq.gz  +  *_R2_val_2.fq.gz
#
# Output (written IN PLACE so step 2 picks them up unchanged)
#   Same filenames; originals are overwritten. Per-sample marker
#     ${TRIM_DIR}/.${sample}_polya.done
#   makes the step idempotent (re-runs skip completed samples).
#   Cutadapt logs: ${TRIM_DIR}/${sample}_polya.log
#   FastQC on the new files is refreshed in ${TRIM_DIR}/fastqc/.
#
# Cutadapt strategy
#   -a "A{100}"  → trims any A-run reaching the 3′ end of R1/R2
#   -g "T{100}"  → trims any T-run starting from the 5′ end of R1/R2
#   -O $POLYA_MIN_STRETCH  effective minimum run length before trimming fires
#   -e $POLYA_ERROR_RATE   mismatch tolerance
#   --minimum-length $TRIM_MIN_LENGTH   discard reads too short to map
#
# Toggle: RUN_POLYA_TRIM in 0_config.sh (default 0; set 1 per project)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "1.1_polya_trim"

check_cmd cutadapt
check_cmd fastqc

POLYA_MIN_STRETCH="${POLYA_MIN_STRETCH:-10}"
POLYA_ERROR_RATE="${POLYA_ERROR_RATE:-0.1}"

echo "=== STEP 1.1: poly-A / poly-T trimming (cutadapt) ==="
echo "[INFO] TRIM_DIR:            $TRIM_DIR"
echo "[INFO] Mode:                $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"
echo "[INFO] POLYA_MIN_STRETCH:   $POLYA_MIN_STRETCH (effective min A/T run length)"
echo "[INFO] POLYA_ERROR_RATE:    $POLYA_ERROR_RATE"
echo "[INFO] --minimum-length:    $TRIM_MIN_LENGTH"

if [[ ! -d "$TRIM_DIR" ]]; then
  echo "[ERROR] TRIM_DIR not found: $TRIM_DIR (run step 1 first)"
  exit 1
fi

if [[ "${SE:-0}" -eq 1 ]]; then
  inputs=( "${TRIM_DIR}"/*_trimmed.fq.gz )
  if (( ${#inputs[@]} == 0 )); then
    # Trimmed files are deleted by step 2 when DELETE_TRIMMED_AFTER_ALIGN=1.
    # If every expected sample already has a .polya.done marker, this is a re-run no-op.
    shopt -s nullglob; markers=( "${TRIM_DIR}"/.*_polya.done ); shopt -u nullglob
    if (( ${#markers[@]} > 0 )); then
      echo "[INFO] No *_trimmed.fq.gz in $TRIM_DIR, but ${#markers[@]} polya.done marker(s) present — already trimmed. Nothing to do."
      exit 0
    fi
    echo "[ERROR] No *_trimmed.fq.gz in $TRIM_DIR (SE mode; run step 1 first)"
    exit 1
  fi
  for R1 in "${inputs[@]}"; do
    base="$(basename "$R1" .fq.gz)"            # <sample>_trimmed
    sample="${base%_trimmed}"
    marker="${TRIM_DIR}/.${sample}_polya.done"
    cutlog="${TRIM_DIR}/${sample}_polya.log"

    if [[ -f "$marker" ]]; then
      echo "[STEP1.1] $sample — marker present, skipping"
      continue
    fi

    tmp_out="${R1}.polya.tmp.fq.gz"
    echo "[STEP1.1] $sample (SE) → cutadapt poly-A/poly-T"
    cutadapt \
      -a "A{100}" \
      -g "T{100}" \
      -O "$POLYA_MIN_STRETCH" \
      -e "$POLYA_ERROR_RATE" \
      --minimum-length "$TRIM_MIN_LENGTH" \
      --cores "$TRIM_CPU" \
      -o "$tmp_out" \
      "$R1" > "$cutlog"

    mv "$tmp_out" "$R1"
    : > "$marker"
  done
else
  inputs=( "${TRIM_DIR}"/*_R1_val_1.fq.gz )
  if (( ${#inputs[@]} == 0 )); then
    shopt -s nullglob; markers=( "${TRIM_DIR}"/.*_polya.done ); shopt -u nullglob
    if (( ${#markers[@]} > 0 )); then
      echo "[INFO] No *_R1_val_1.fq.gz in $TRIM_DIR, but ${#markers[@]} polya.done marker(s) present — already trimmed. Nothing to do."
      exit 0
    fi
    echo "[ERROR] No *_R1_val_1.fq.gz in $TRIM_DIR (PE mode; run step 1 first)"
    exit 1
  fi
  for R1 in "${inputs[@]}"; do
    R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing mate for $R1 (expected $R2). Skipping."
      continue
    fi
    sample="$(basename "$R1" _R1_val_1.fq.gz)"
    marker="${TRIM_DIR}/.${sample}_polya.done"
    cutlog="${TRIM_DIR}/${sample}_polya.log"

    if [[ -f "$marker" ]]; then
      echo "[STEP1.1] $sample — marker present, skipping"
      continue
    fi

    tmp_r1="${R1}.polya.tmp.fq.gz"
    tmp_r2="${R2}.polya.tmp.fq.gz"
    echo "[STEP1.1] $sample (PE) → cutadapt poly-A/poly-T (R1+R2)"
    cutadapt \
      -a "A{100}" -g "T{100}" \
      -A "A{100}" -G "T{100}" \
      -O "$POLYA_MIN_STRETCH" \
      -e "$POLYA_ERROR_RATE" \
      --minimum-length "$TRIM_MIN_LENGTH" \
      --cores "$TRIM_CPU" \
      -o "$tmp_r1" -p "$tmp_r2" \
      "$R1" "$R2" > "$cutlog"

    mv "$tmp_r1" "$R1"
    mv "$tmp_r2" "$R2"
    : > "$marker"
  done
fi

###############################################################################
# Refresh FastQC on the post-poly-A outputs so step-5 MultiQC sees the final state
###############################################################################

if [[ "${SE:-0}" -eq 1 ]]; then
  trimmed=( "${TRIM_DIR}"/*_trimmed.fq.gz )
else
  trimmed=( "${TRIM_DIR}"/*_val_*.fq.gz )
fi

fastqc_dir="${TRIM_DIR}/fastqc"
mkdir -p "$fastqc_dir"
if (( ${#trimmed[@]} > 0 )); then
  rm -f "${fastqc_dir}"/*.html "${fastqc_dir}"/*.zip
  echo "[STEP1.1] FastQC on poly-A-trimmed reads → $fastqc_dir"
  fastqc --outdir "$fastqc_dir" --threads "$TRIM_CPU" "${trimmed[@]}"
fi

echo "=== STEP 1.1 complete ==="
