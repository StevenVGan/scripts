#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.2_peak_homer.sh — HOMER findPeaks for PRO-seq (default -style groseq)
#
# Purpose
#   HOMER-based feature calling + annotation. For PRO-seq the default style is
#   "groseq", which detects nascent transcript intervals from stranded tag
#   directories — not the TF/histone peak concept used for ChIP-style data.
#
#   The input groups file format is the same as csRNA/cutrun (for compatibility
#   and reuse of helpers), but:
#     - controls are almost always absent for PRO-seq → rows with "-"/none/NA
#       in column 2 work as expected;
#     - the "type" column in column 4 is IGNORED for PRO-seq — HOMER_STYLE from
#       0_config.sh (default "groseq") is used instead.
#
#   If PEAKCALL_GROUPS_FILE does not exist, this script auto-creates a default
#   one (one row per _sorted.bam, no control, type=TF) so you can run without
#   editing anything.
#
# Flow per group
#   1) findPeaks $tag_dir -style $HOMER_STYLE [-i $ctrl_tag_dir] -o $out_file
#   2) annotatePeaks.pl → ${name}.annotatePeaks.txt
#   3) Save findPeaks header lines to ${name}_${suffix}_stats.txt; remove the
#      full peak .txt (findPeaks output is large and not used after annotation).
#   4) Filter annotatePeaks by blacklist.
#
# Outputs
#   - ${HOMER_PEAK_DIR}/${name}_${suffix}_stats.txt  (headers only, for QC)
#   - ${HOMER_PEAK_DIR}/${name}.annotatePeaks.txt    (annotated, blacklist-filtered)
#
#   suffix by style:  groseq → "transcripts", histone → "regions", else → "peaks"
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.2_peak_homer"

check_cmd findPeaks
check_cmd annotatePeaks.pl
check_cmd bedtools

echo "=== STEP 4.2: HOMER findPeaks + annotatePeaks + blacklist filter ==="
echo "[INFO] HOMER_STYLE:      ${HOMER_STYLE:-groseq}"
echo "[INFO] Groups file:      $PEAKCALL_GROUPS_FILE"
echo "[INFO] Tag dir:          $TAG_DIR"
echo "[INFO] HOMER peak dir:   $HOMER_PEAK_DIR"
echo "[INFO] GENOME:           $GENOME"
echo "[INFO] BLACKLIST:        $BLACKLIST"

# Auto-create groups file if missing (PRO-seq default: one row per BAM, no control).
if [[ ! -f "$PEAKCALL_GROUPS_FILE" ]]; then
  echo "[INFO] Groups file not found — auto-generating: $PEAKCALL_GROUPS_FILE"
  : > "$PEAKCALL_GROUPS_FILE"
  for bam in "${BAM_DIR}"/*_sorted.bam; do
    [[ -f "$bam" ]] || continue
    b="$(basename "$bam")"
    n="$(basename "$bam" _sorted.bam)"
    printf '%s\t-\t%s\tTF\n' "$b" "$n" >> "$PEAKCALL_GROUPS_FILE"
  done
fi

if [[ ! -s "$PEAKCALL_GROUPS_FILE" ]]; then
  echo "[ERROR] Groups file is empty: $PEAKCALL_GROUPS_FILE"
  exit 1
fi

if [[ ! -f "$BLACKLIST" ]]; then
  echo "[WARN] Blacklist not found; will skip blacklist filtering on annotatePeaks."
  SKIP_BLACKLIST=1
else
  SKIP_BLACKLIST=0
fi

mkdir -p "$HOMER_PEAK_DIR"

# File naming suffix by findPeaks style.
style_suffix() {
  case "$1" in
    groseq)  echo "transcripts" ;;
    histone) echo "regions" ;;
    *)       echo "peaks" ;;
  esac
}

while IFS=$'\t' read -r ip_spec ctrl_spec name type _rest; do
  [[ -z "${ip_spec:-}" ]] && continue
  [[ "$ip_spec" =~ ^# ]] && continue

  if [[ -z "${name:-}" ]]; then name="$(basename "$ip_spec" .bam)"; fi

  # Choose style: HOMER_STYLE from config always wins for PRO-seq, unless the
  # row's type column explicitly says Histone (escape hatch for mixed analyses).
  style="${HOMER_STYLE:-groseq}"
  type_lower="$(echo "${type:-}" | tr '[:upper:]' '[:lower:]')"
  [[ "$type_lower" == "histone" ]] && style="histone"

  peak_suffix="$(style_suffix "$style")"

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

  out_file="${HOMER_PEAK_DIR}/${name}_${peak_suffix}.txt"
  ann_file="${HOMER_PEAK_DIR}/${name}.annotatePeaks.txt"
  stats_file="${HOMER_PEAK_DIR}/${name}_${peak_suffix}_stats.txt"

  echo "[STEP4.2] $name  style=$style  ctrl=$([ $use_control -eq 1 ] && echo yes || echo no)"

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

  # 4) Filter annotatePeaks by blacklist (col1=PeakID, col2=Chr, col3=Start, col4=End)
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
