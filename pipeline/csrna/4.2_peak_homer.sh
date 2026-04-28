#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.2_peak_homer.sh
#
# Purpose
#   HOMER peak calling (findPeaks) + annotation (annotatePeaks) in one step.
#   Uses peakcall_groups: TF → factor style, Histone → histone style.
#   Keeps only findPeaks headers as a small QC stats file; removes full peak .txt.
#   Final output: annotatePeaks table, filtered by blacklist.
#
# Flow per group
#   1) findPeaks → *_peaks.txt or *_regions.txt (by type)
#   2) annotatePeaks.pl → ${name}.annotatePeaks.txt
#   3) Save findPeaks header lines to ${name}_peaks_stats.txt or _regions_stats.txt
#   4) Remove the full *_peaks.txt / *_regions.txt
#   5) Filter .annotatePeaks.txt by blacklist (remove overlapping rows)
#
# Outputs
#   - ${HOMER_PEAK_DIR}/${name}_peaks_stats.txt or _regions_stats.txt (headers only, for QC)
#   - ${HOMER_PEAK_DIR}/${name}.annotatePeaks.txt (annotated peaks, blacklist-filtered)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.2_peak_homer"

check_cmd findPeaks
check_cmd annotatePeaks.pl
check_cmd bedtools

echo "=== STEP 4.2: HOMER findPeaks + annotatePeaks + blacklist filter ==="
echo "[INFO] Groups file:   $PEAKCALL_GROUPS_FILE"
echo "[INFO] Tag dir:       $TAG_DIR"
echo "[INFO] HOMER peak dir: $HOMER_PEAK_DIR"
echo "[INFO] GENOME:        $GENOME"
echo "[INFO] BLACKLIST:     $BLACKLIST"

if [[ ! -f "$PEAKCALL_GROUPS_FILE" ]]; then
  echo "[ERROR] Groups file not found: $PEAKCALL_GROUPS_FILE"
  exit 1
fi

if [[ ! -f "$BLACKLIST" ]]; then
  echo "[WARN] Blacklist not found; will skip blacklist filtering on annotatePeaks."
  SKIP_BLACKLIST=1
else
  SKIP_BLACKLIST=0
fi

mkdir -p "$HOMER_PEAK_DIR"

while IFS=$'\t' read -r ip_spec ctrl_spec name type _rest; do
  [[ -z "${ip_spec:-}" ]] && continue
  [[ "$ip_spec" =~ ^# ]] && continue

  if [[ -z "${name:-}" ]]; then name="$(basename "$ip_spec" .bam)"; fi
  if [[ -z "${type:-}" ]]; then type="TF"; fi
  type_lower="$(echo "$type" | tr '[:upper:]' '[:lower:]')"

  if [[ -f "$ip_spec" ]]; then
    ip_file="$ip_spec"
  elif [[ -f "${BAM_DIR}/${ip_spec}" ]]; then
    ip_file="${BAM_DIR}/${ip_spec}"
  else
    echo "[WARN] IP BAM not found for '$ip_spec'. Skipping $name"
    continue
  fi

  ip_base="$(basename "$ip_file" _sorted.bam)"
  tag_dir="${TAG_DIR}/${ip_base}"
  if [[ ! -d "$tag_dir" ]]; then
    echo "[WARN] Tag dir not found: $tag_dir. Skipping $name"
    continue
  fi

  use_control=1
  if [[ -z "${ctrl_spec:-}" ]] || [[ "$ctrl_spec" == "-" ]] || [[ "$ctrl_spec" =~ ^(none|NONE|NA)$ ]]; then
    use_control=0
  fi

  ctrl_tag_dir=""
  if [[ "$use_control" -eq 1 ]]; then
    if [[ -f "$ctrl_spec" ]]; then
      ctrl_bam="$ctrl_spec"
    elif [[ -f "${BAM_DIR}/${ctrl_spec}" ]]; then
      ctrl_bam="${BAM_DIR}/${ctrl_spec}"
    else
      echo "[WARN] Control not found for '$ctrl_spec'. Skipping $name"
      continue
    fi
    ctrl_base="$(basename "$ctrl_bam" _sorted.bam)"
    ctrl_tag_dir="${TAG_DIR}/${ctrl_base}"
    if [[ ! -d "$ctrl_tag_dir" ]]; then
      echo "[WARN] Control tag dir not found: $ctrl_tag_dir. Skipping $name"
      continue
    fi
  fi

  if [[ "$type_lower" == "histone" ]]; then
    style="histone"
    peak_suffix="regions"
  else
    style="factor"
    peak_suffix="peaks"
  fi
  out_file="${HOMER_PEAK_DIR}/${name}_${peak_suffix}.txt"
  ann_file="${HOMER_PEAK_DIR}/${name}.annotatePeaks.txt"
  stats_file="${HOMER_PEAK_DIR}/${name}_${peak_suffix}_stats.txt"

  echo "[STEP4.2] $name  style=$style  -> annotate -> stats only -> filter annotatePeaks"

  # 1) findPeaks
  if [[ "$use_control" -eq 1 ]]; then
    findPeaks "$tag_dir" -style "$style" -o "$out_file" -i "$ctrl_tag_dir"
  else
    findPeaks "$tag_dir" -style "$style" -o "$out_file"
  fi

  [[ ! -f "$out_file" ]] && continue

  # 2) annotatePeaks
  annotatePeaks.pl "$out_file" "$GENOME" > "$ann_file"

  # 3) Keep only findPeaks headers for QC stats, then remove full peak file
  grep '^#' "$out_file" > "$stats_file"
  rm -f "$out_file"

  # 4) Filter annotatePeaks by blacklist (annotatePeaks: col1=PeakID, col2=Chr, col3=Start, col4=End)
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
