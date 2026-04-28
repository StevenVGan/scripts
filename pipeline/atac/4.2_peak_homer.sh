#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.2_peak_homer.sh — HOMER findPeaks for ATAC (-style dnase, no control)
#
# Per group:
#   1) findPeaks $tag_dir -style dnase -o $out_file
#   2) annotatePeaks.pl $out_file $GENOME > ${name}.annotatePeaks.txt
#   3) Save findPeaks header lines to ${name}_peaks_stats.txt; remove full peak file
#   4) Filter annotatePeaks by blacklist
#
# Differences vs cutrun: -style dnase always (no TF/Histone branch),
# no control tag dir lookup, and tag_dir derivation strips _sorted.shifted.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.2_peak_homer"

check_cmd findPeaks
check_cmd annotatePeaks.pl
check_cmd bedtools

echo "=== STEP 4.2: HOMER findPeaks -style dnase + annotatePeaks + blacklist filter ==="
echo "[INFO] Groups file:    $PEAKCALL_GROUPS_FILE"
echo "[INFO] Tag dir:        $TAG_DIR"
echo "[INFO] HOMER peak dir: $HOMER_PEAK_DIR"
echo "[INFO] GENOME:         $GENOME"
echo "[INFO] BLACKLIST:      $BLACKLIST"

if [[ ! -f "$PEAKCALL_GROUPS_FILE" ]]; then
  echo "[ERROR] Groups file not found: $PEAKCALL_GROUPS_FILE"
  exit 1
fi

if [[ ! -f "$BLACKLIST" ]]; then
  echo "[WARN] Blacklist not found; will skip blacklist filtering."
  SKIP_BLACKLIST=1
else
  SKIP_BLACKLIST=0
fi

mkdir -p "$HOMER_PEAK_DIR"

while IFS=$'\t' read -r ip_spec ctrl_spec name type _rest; do
  [[ -z "${ip_spec:-}" ]] && continue
  [[ "$ip_spec" =~ ^# ]] && continue

  if [[ -z "${name:-}" ]]; then name="$(basename "$ip_spec" .bam)"; fi

  if [[ -f "$ip_spec" ]]; then
    ip_file="$ip_spec"
  elif [[ -f "${BAM_DIR}/${ip_spec}" ]]; then
    ip_file="${BAM_DIR}/${ip_spec}"
  else
    echo "[WARN] IP BAM not found for '$ip_spec'. Skipping $name"
    continue
  fi

  # Strip _sorted.shifted.bam (preferred) or _sorted.bam (fallback) to get tag-dir basename
  ip_base="$(basename "$ip_file")"
  ip_base="${ip_base%_sorted.shifted.bam}"
  ip_base="${ip_base%_sorted.bam}"
  ip_base="${ip_base%.bam}"
  tag_dir="${TAG_DIR}/${ip_base}"
  if [[ ! -d "$tag_dir" ]]; then
    echo "[WARN] Tag dir not found: $tag_dir. Skipping $name"
    continue
  fi

  if [[ -n "${ctrl_spec:-}" && "$ctrl_spec" != "-" && ! "$ctrl_spec" =~ ^(none|NONE|NA)$ ]]; then
    echo "[WARN] $name: control '$ctrl_spec' specified but ATAC pipeline ignores controls."
  fi

  out_file="${HOMER_PEAK_DIR}/${name}_peaks.txt"
  ann_file="${HOMER_PEAK_DIR}/${name}.annotatePeaks.txt"
  stats_file="${HOMER_PEAK_DIR}/${name}_peaks_stats.txt"

  echo "[STEP4.2] $name  style=dnase  tag_dir=$tag_dir"

  # 1) findPeaks (no control)
  findPeaks "$tag_dir" -style dnase -o "$out_file"
  [[ ! -f "$out_file" ]] && continue

  # 2) annotatePeaks
  annotatePeaks.pl "$out_file" "$GENOME" > "$ann_file"

  # 3) Headers-only stats; drop full peak file
  grep '^#' "$out_file" > "$stats_file"
  rm -f "$out_file"

  # 4) Blacklist filter on annotatePeaks (col1=PeakID, col2=Chr, col3=Start, col4=End)
  if [[ "$SKIP_BLACKLIST" -eq 1 ]] || [[ ! -f "$ann_file" ]]; then
    continue
  fi
  kept_tmp=$(mktemp)
  ann_tmp=$(mktemp)
  tail -n +2 "$ann_file" | awk -v OFS='\t' 'NF>=4 {print $2,$3,$4,$1}' | sort -k1,1 -k2,2n \
    | bedtools intersect -a stdin -b "$BLACKLIST" -v | cut -f4 > "$kept_tmp"
  head -1 "$ann_file" > "$ann_tmp"
  tail -n +2 "$ann_file" | awk -v kept="$kept_tmp" 'BEGIN{while((getline < kept)>0) k[$1]=1} $1 in k{print}' >> "$ann_tmp"
  mv "$ann_tmp" "$ann_file"
  rm -f "$kept_tmp"
done < "$PEAKCALL_GROUPS_FILE"

echo "=== STEP 4.2 (HOMER) complete ==="
