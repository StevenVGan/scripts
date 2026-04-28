#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 3_homer_tags.sh
#
# Purpose
#   Step 3 of the pipeline: create HOMER tag directories from aligned BAMs.
#   Tag directories are used by HOMER tools for downstream analysis/visualization.
#
# Inputs
#   - Sorted BAMs in ${BAM_DIR}:
#       *_sorted.bam
#
# Outputs
#   - HOMER tag directories in ${TAG_DIR}:
#       ${TAG_DIR}/${sample_name}/  (one directory per sample)
#
# Behavior
#   - Iterates over all *_sorted.bam files.
#   - Skips samples if their tag directory already exists.
#   - Uses makeTagDirectory with "-tbp 1" (format auto-detected; HOMER converts BAM on the fly).
#
# Logging
#   - Writes a timestamped log to ${LOG_DIR}/3_homer_tags_*.log
#
# Requirements
#   - makeTagDirectory (HOMER) on PATH.
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "3_homer_tags"

check_cmd makeTagDirectory

echo "=== STEP 3: HOMER makeTagDirectory from sorted BAMs ==="
echo "[INFO] BAM_DIR: $BAM_DIR"
echo "[INFO] TAG_DIR: $TAG_DIR"

sorted_bams=( "${BAM_DIR}"/*_sorted.bam )
if (( ${#sorted_bams[@]} == 0 )); then
  echo "[ERROR] No *_sorted.bam found in $BAM_DIR"
  exit 1
fi

for bam in "${sorted_bams[@]}"; do
  sample_name="$(basename "$bam" _sorted.bam)"
  outdir="${TAG_DIR}/${sample_name}"

  if [[ -d "$outdir" ]]; then
    echo "[STEP3] Tag dir exists for $sample_name, skipping: $outdir"
    continue
  fi

  echo "[STEP3] Building tag dir: $sample_name"
  makeTagDirectory "$outdir" "$bam" -tbp 1
done

echo "=== STEP 3 complete ==="
