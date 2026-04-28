#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 3_homer_tags.sh — HOMER tag directories from Tn5-shifted BAMs
#
# Difference vs cutrun: globs *_sorted.shifted.bam (not *_sorted.bam) and
# strips the right suffix so tag dirs are named ${sample}/ to match the
# basename used in peakcall_groups.tsv.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "3_homer_tags"

check_cmd makeTagDirectory

echo "=== STEP 3: HOMER makeTagDirectory from Tn5-shifted BAMs ==="
echo "[INFO] BAM_DIR: $BAM_DIR"
echo "[INFO] TAG_DIR: $TAG_DIR"

shifted_bams=( "${BAM_DIR}"/*_sorted.shifted.bam )
if (( ${#shifted_bams[@]} == 0 )); then
  echo "[ERROR] No *_sorted.shifted.bam found in $BAM_DIR (run 2_bowtie2 first)"
  exit 1
fi

for bam in "${shifted_bams[@]}"; do
  sample_name="$(basename "$bam" _sorted.shifted.bam)"
  outdir="${TAG_DIR}/${sample_name}"

  if [[ -d "$outdir" ]]; then
    echo "[STEP3] Tag dir exists for $sample_name, skipping: $outdir"
    continue
  fi

  echo "[STEP3] Building tag dir: $sample_name"
  makeTagDirectory "$outdir" "$bam" -tbp 1
done

echo "=== STEP 3 complete ==="
