#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_all.sh — ATAC-seq upstream pipeline orchestrator
#
#   1)   1_trim_qc.sh      Trim Galore (--nextera) + FastQC
#   2)   2_bowtie2.sh      bowtie2 (-X 2000) + dedup + filter + Tn5 shift + bigWig
#   3)   3_homer_tags.sh   makeTagDirectory from *_sorted.shifted.bam
#   4.1) 4.1_peak_macs3.sh MACS3 callpeak (--nomodel --shift -75 --extsize 150)
#   4.2) 4.2_peak_homer.sh findPeaks -style dnase
#   5)   5_qc.sh           preseq + deepTools + fragment-size + TSS enrichment + MultiQC
#
# Toggles (env): RUN_TRIM, RUN_BOWTIE2, RUN_HOMER_TAGS, RUN_PEAK_MACS3,
#                RUN_PEAK_HOMER, RUN_QC.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

echo "=== run_all.sh starting (ATAC) ==="
echo "[INFO] BASE: $BASE"
echo "[INFO] Toggles: TRIM=$RUN_TRIM BOWTIE2=$RUN_BOWTIE2 TAGS=$RUN_HOMER_TAGS MACS3=$RUN_PEAK_MACS3 HOMER_PEAK=$RUN_PEAK_HOMER QC=$RUN_QC"

if [[ "$RUN_TRIM" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/1_trim_qc.sh"
else
  echo "=== Skipping STEP 1 (RUN_TRIM=0) ==="
fi

if [[ "$RUN_BOWTIE2" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/2_bowtie2.sh"
else
  echo "=== Skipping STEP 2 (RUN_BOWTIE2=0) ==="
fi

if [[ "$RUN_HOMER_TAGS" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/3_homer_tags.sh"
else
  echo "=== Skipping STEP 3 (RUN_HOMER_TAGS=0) ==="
fi

if [[ "$RUN_PEAK_MACS3" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/4.1_peak_macs3.sh"
else
  echo "=== Skipping STEP 4.1 (RUN_PEAK_MACS3=0) ==="
fi

if [[ "$RUN_PEAK_HOMER" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/4.2_peak_homer.sh"
else
  echo "=== Skipping STEP 4.2 (RUN_PEAK_HOMER=0) ==="
fi

if [[ "$RUN_QC" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/5_qc.sh"
else
  echo "=== Skipping STEP 5 (RUN_QC=0) ==="
fi

echo "=== run_all.sh finished successfully ==="
