#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 3_homer_tags.sh — HOMER tag directories for PRO-seq (strand-flip aware)
#
# Forked from csRNA/3_homer_tags.sh. Adds -flip to makeTagDirectory when
# PROSEQ_FLIP_STRAND=1 so that all downstream HOMER tools (findPeaks -style
# groseq, analyzeRNA, annotatePeaks) see the nascent RNA strand rather than
# the read strand.
#
# Inputs
#   - Sorted BAMs in ${BAM_DIR}: *_sorted.bam
#
# Outputs
#   - ${TAG_DIR}/${sample_name}/  (one tag directory per sample)
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "3_homer_tags"

check_cmd makeTagDirectory

echo "=== STEP 3: HOMER makeTagDirectory from sorted BAMs ==="
echo "[INFO] BAM_DIR:            $BAM_DIR"
echo "[INFO] TAG_DIR:            $TAG_DIR"
echo "[INFO] HOMER_SS_PE:        ${HOMER_SS_PE:-0}"
echo "[INFO] PROSEQ_FLIP_STRAND: ${PROSEQ_FLIP_STRAND:-0} (adds -flip when 1)"

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
  homer_extra=( -tbp 1 )
  [[ "${HOMER_SS_PE:-0}"         -eq 1 ]] && homer_extra+=( -sspe )
  [[ "${PROSEQ_FLIP_STRAND:-0}"  -eq 1 ]] && homer_extra+=( -flip )
  makeTagDirectory "$outdir" "$bam" "${homer_extra[@]}"
done

echo "=== STEP 3 complete ==="
