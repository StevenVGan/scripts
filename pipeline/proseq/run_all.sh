#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_all.sh — PRO-seq pipeline orchestrator
#
# Steps (fork of csRNA with PRO-seq-specific strand flip + pausing module):
#   1)   Trim Galore (auto R1 adapter, --nextseq polyG, --adapter2 SmallRNA-5p
#        for R2 readthrough) + FastQC
#   1.1) Optional poly-A / poly-T trim (cutadapt; enable with RUN_POLYA_TRIM=1
#        per project when FastQC flags poly-A content, e.g. mRNA carryover)
#   2)   bowtie2 --very-sensitive-local → sorted BAM + stranded bigWigs
#        (bigWig labels swapped when PROSEQ_FLIP_STRAND=1 so _fwd.bw = RNA on +)
#   3)   HOMER makeTagDirectory (-sspe for PE, -flip for PRO-seq strand)
#   4.1) MACS3 peak calling (off by default; rarely appropriate for nascent RNA)
#   4.2) HOMER findPeaks -style groseq → annotated blacklist-filtered transcripts
#   4.3) Pausing index + divergent eRNA pairs
#   5)   MultiQC + QC (preseq, plotFingerprint, PCA/Correlation, pausing summary)
#
# Toggles: RUN_TRIM, RUN_POLYA_TRIM, RUN_BOWTIE2, RUN_HOMER_TAGS, RUN_PEAK_MACS3,
#          RUN_PEAK_HOMER, RUN_PAUSING, RUN_QC
#
# Usage
#   ./run_all.sh                                        # full pipeline
#   RUN_TRIM=0 RUN_BOWTIE2=0 RUN_PEAK_HOMER=1 ./run_all.sh   # resume from peak calling
#   SE=1 ./run_all.sh                                   # single-end project
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

echo "=== run_all.sh (PRO-seq) starting ==="
echo "[INFO] BASE: $BASE"
echo "[INFO] Toggles: TRIM=$RUN_TRIM POLYA=$RUN_POLYA_TRIM BOWTIE2=$RUN_BOWTIE2 TAGS=$RUN_HOMER_TAGS MACS3=$RUN_PEAK_MACS3 HOMER_PEAK=$RUN_PEAK_HOMER PAUSING=$RUN_PAUSING QC=$RUN_QC"
echo "[INFO] Mode:   SE=$SE  PROSEQ_FLIP_STRAND=$PROSEQ_FLIP_STRAND  HOMER_STYLE=$HOMER_STYLE"

if [[ "$RUN_TRIM" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/1_trim_qc.sh"
else
  echo "=== Skipping STEP 1 (RUN_TRIM=0) ==="
fi

if [[ "$RUN_POLYA_TRIM" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/1.1_polya_trim.sh"
else
  echo "=== Skipping STEP 1.1 (RUN_POLYA_TRIM=0) ==="
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

if [[ "$RUN_PAUSING" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/4.3_pausing_divergent.sh"
else
  echo "=== Skipping STEP 4.3 (RUN_PAUSING=0) ==="
fi

if [[ "$RUN_QC" -eq 1 ]]; then
  bash "${SCRIPT_DIR}/5_qc.sh"
else
  echo "=== Skipping STEP 5 (RUN_QC=0) ==="
fi

echo "=== run_all.sh (PRO-seq) finished successfully ==="
