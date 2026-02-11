#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_all.sh
#
# Purpose
#   Orchestrator script that runs the pipeline steps in order:
#     1) trimming + FastQC
#     2) bowtie2 alignment + BAM/track generation
#     3) HOMER tag directories
#     4.1) MACS3 peak calling (4.1_peak_macs3.sh)
#     4.2) HOMER peak calling + annotation + blacklist filter (4.2_peak_homer.sh)
#     5) QC + MultiQC (5_qc.sh: preseq, run_spp.R, plotFingerprint, PCA/Correlation, optional MACS3 summary)
#
# Configuration
#   - Sources 0_config.sh. Toggles: RUN_TRIM, RUN_BOWTIE2, RUN_HOMER_TAGS, RUN_PEAK_MACS3, RUN_PEAK_HOMER, RUN_QC.
#
# Usage
#   - Run everything: ./run_all.sh
#   - Run only peak calling: RUN_TRIM=0 RUN_BOWTIE2=0 RUN_HOMER_TAGS=0 RUN_PEAK_MACS3=1 RUN_PEAK_HOMER=1 ./run_all.sh
#   - Run only QC: RUN_TRIM=0 RUN_BOWTIE2=0 RUN_HOMER_TAGS=0 RUN_PEAK_MACS3=0 RUN_PEAK_HOMER=0 RUN_QC=1 ./run_all.sh
#
# Behavior
#   - Stops immediately if a called step exits with a non-zero status
#     (because each step uses set -euo pipefail).
#
# Logging
#   - Each step script writes its own timestamped log in ${LOG_DIR}.
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

echo "=== run_all.sh starting ==="
echo "[INFO] BASE: $BASE"

# You can override toggles via env vars, e.g.:
# RUN_BOWTIE2=0 ./run_all.sh
RUN_TRIM="${RUN_TRIM:-1}"
RUN_BOWTIE2="${RUN_BOWTIE2:-1}"
RUN_HOMER_TAGS="${RUN_HOMER_TAGS:-1}"
RUN_PEAK_MACS3="${RUN_PEAK_MACS3:-1}"
RUN_PEAK_HOMER="${RUN_PEAK_HOMER:-1}"
RUN_QC="${RUN_QC:-1}"

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
