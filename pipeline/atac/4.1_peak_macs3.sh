#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.1_peak_macs3.sh — MACS3 peak calling for ATAC (cut-site mode)
#
# ATAC-specific MACS3 invocation (forced for SE and PE):
#   macs3 callpeak -t $ip_bam -f $MACS3_FORMAT -g hs -q $MACS3_FDR \
#       --nomodel --shift $MACS3_SHIFT --extsize $MACS3_EXTSIZE \
#       --keep-dup $MACS3_KEEP_DUP --outdir $MACS3_DIR --name $name
#
# - No --broad branch (ATAC peaks are point-source).
# - No control branch (ATAC has no input/IgG).
# - Treats every read as a Tn5 cut site (works on shifted PE too — MACS3 ignores
#   pair info when format=BAM, which is what we want for cut-site peak calling).
#
# peakcall_groups.tsv columns: ip_bam <TAB> control_bam <TAB> name <TAB> type
#   - control_bam should be "-" / "none" / "NA" — non-blank values are ignored with a warning.
#   - type column is informational; this script always calls narrow.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.1_peak_macs3"

check_cmd_string "$MACS3_CMD"
check_cmd bedtools
check_cmd annotatePeaks.pl

echo "=== STEP 4.1: MACS3 (ATAC cut-site mode) + blacklist filter + annotatePeaks ==="
echo "[INFO] Groups file: $PEAKCALL_GROUPS_FILE"
echo "[INFO] MACS3 outdir: $MACS3_DIR"
echo "[INFO] BLACKLIST:    $BLACKLIST"
echo "[INFO] Params: --nomodel --shift $MACS3_SHIFT --extsize $MACS3_EXTSIZE --keep-dup $MACS3_KEEP_DUP --format $MACS3_FORMAT"

if [[ ! -f "$PEAKCALL_GROUPS_FILE" ]]; then
  echo "[ERROR] Groups file not found: $PEAKCALL_GROUPS_FILE"
  exit 1
fi

if [[ ! -f "$BLACKLIST" ]]; then
  echo "[WARN] Blacklist not found: $BLACKLIST; will skip blacklist filtering."
  SKIP_BLACKLIST=1
else
  SKIP_BLACKLIST=0
fi

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

  if [[ -n "${ctrl_spec:-}" && "$ctrl_spec" != "-" && ! "$ctrl_spec" =~ ^(none|NONE|NA)$ ]]; then
    echo "[WARN] $name: control '$ctrl_spec' specified but ATAC pipeline ignores controls."
  fi

  echo "[STEP4.1] $name  IP: $ip_file"

  # shellcheck disable=SC2206
  macs_cmd=( $MACS3_CMD callpeak
    -t "$ip_file"
    -g "$MACS3_GENOMESIZE"
    -q "$MACS3_FDR"
    --format "$MACS3_FORMAT"
    --nomodel
    --shift "$MACS3_SHIFT"
    --extsize "$MACS3_EXTSIZE"
    --keep-dup "$MACS3_KEEP_DUP"
    --outdir "$MACS3_DIR"
    --name "$name"
  )
  "${macs_cmd[@]}"

  narrowPeak_file="${MACS3_DIR}/${name}_peaks.narrowPeak"
  filtered_bed="${MACS3_DIR}/${name}_filtered.bed"

  if [[ "$SKIP_BLACKLIST" -eq 0 && -f "$narrowPeak_file" ]]; then
    bedtools intersect -a "$narrowPeak_file" -b "$BLACKLIST" -v > "$filtered_bed"
  fi

  bed_to_annotate=""
  if [[ -f "$filtered_bed" ]]; then
    bed_to_annotate="$filtered_bed"
  elif [[ -f "$narrowPeak_file" ]]; then
    bed_to_annotate="$narrowPeak_file"
  fi
  if [[ -n "$bed_to_annotate" ]]; then
    out_ann="${MACS3_DIR}/${name}.annotatePeaks.txt"
    echo "[STEP4.1] Annotate: $bed_to_annotate -> $out_ann"
    annotatePeaks.pl "$bed_to_annotate" "$GENOME" > "$out_ann"
  else
    echo "[WARN] No peak file found for $name."
  fi
done < "$PEAKCALL_GROUPS_FILE"

echo "=== STEP 4.1 (MACS3) complete ==="
