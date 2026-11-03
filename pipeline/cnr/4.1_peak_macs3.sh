#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.1_peak_macs3.sh
#
# Purpose
#   Step 4.1: peak calling with MACS3 and (optional) blacklist filtering.
#   Uses shared peak-call groups file (see PEAKCALL_GROUPS_FILE).
#
# Inputs
#   - IP/control BAMs described by ${PEAKCALL_GROUPS_FILE}, format:
#       ip_bam   control_bam   name   type
#     * control_bam can be "-", "none", or "NA" for no control
#     * type: TF or Histone (Histone → broad calling)
#
# Outputs
#   - ${MACS3_DIR}/${name}_peaks.narrowPeak or .broadPeak
#   - ${MACS3_DIR}/${name}_filtered.bed (if PEAK_BLACKLIST exists; default = BLACKLIST)
#   - ${MACS3_DIR}/${name}.annotatePeaks.txt (HOMER annotatePeaks on filtered BED or peak file)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.1_peak_macs3"

# Peak-filter BED (2026-10-26, opt-in): a project may set PEAK_BLACKLIST in its 0_config.sh to a superset of BLACKLIST
# (e.g. v2 + whole chrM/unplaced contigs, CnR_261006 review R-02) that is applied to the peak calls only; bamCoverage and
# the deepTools QC keep BLACKLIST. Unset = BLACKLIST, so every other project is unchanged.
PEAK_BLACKLIST="${PEAK_BLACKLIST:-${BLACKLIST}}"

check_cmd_string "$MACS3_CMD"
check_cmd bedtools
check_cmd annotatePeaks.pl

echo "=== STEP 4.1: MACS3 peak calling + blacklist filtering + annotatePeaks ==="
echo "[INFO] Groups file: $PEAKCALL_GROUPS_FILE"
echo "[INFO] MACS3 outdir: $MACS3_DIR"
echo "[INFO] BLACKLIST:    $BLACKLIST"
echo "[INFO] PEAK_BLACKLIST (peak filter): $PEAK_BLACKLIST"

if [[ ! -f "$PEAKCALL_GROUPS_FILE" ]]; then
  echo "[ERROR] Groups file not found: $PEAKCALL_GROUPS_FILE"
  exit 1
fi

if [[ ! -f "$PEAK_BLACKLIST" && "$PEAK_BLACKLIST" != "$BLACKLIST" ]]; then
  # a project that opted in must not silently fall back to unfiltered calls (e.g. a node without that BED)
  echo "[ERROR] PEAK_BLACKLIST is set but not found: $PEAK_BLACKLIST (build it, or unset PEAK_BLACKLIST to filter with BLACKLIST)"
  exit 1
elif [[ ! -f "$PEAK_BLACKLIST" ]]; then
  echo "[WARN] Peak blacklist not found: $PEAK_BLACKLIST; will skip blacklist filtering."
  SKIP_BLACKLIST=1
else
  SKIP_BLACKLIST=0
fi

while IFS=$'\t' read -r ip_spec ctrl_spec name type _rest; do
  [[ -z "${ip_spec:-}" ]] && continue
  [[ "$ip_spec" =~ ^# ]] && continue

  if [[ -z "${name:-}" ]]; then name="$(basename "$ip_spec" .bam)"; fi
  if [[ -z "${type:-}" ]]; then type="TF"; fi
  type_upper="$(echo "$type" | tr '[:lower:]' '[:upper:]')"

  if [[ -f "$ip_spec" ]]; then
    ip_file="$ip_spec"
  elif [[ -f "${BAM_DIR}/${ip_spec}" ]]; then
    ip_file="${BAM_DIR}/${ip_spec}"
  else
    echo "[WARN] IP BAM not found for '$ip_spec'. Skipping $name"
    continue
  fi

  use_control=1
  if [[ -z "${ctrl_spec:-}" ]] || [[ "$ctrl_spec" == "-" ]] || [[ "$ctrl_spec" =~ ^(none|NONE|NA)$ ]]; then
    use_control=0
  fi

  control_file=""
  if [[ "$use_control" -eq 1 ]]; then
    if [[ -f "$ctrl_spec" ]]; then
      control_file="$ctrl_spec"
    elif [[ -f "${BAM_DIR}/${ctrl_spec}" ]]; then
      control_file="${BAM_DIR}/${ctrl_spec}"
    else
      echo "[WARN] Control BAM not found for '$ctrl_spec'. Skipping $name"
      continue
    fi
  fi

  echo "[STEP4.1] Group: $name  TYPE=$type_upper  IP: $ip_file  CTRL: ${control_file:-(none)}"

  # shellcheck disable=SC2206
  macs_cmd=( $MACS3_CMD callpeak
    -t "$ip_file"
    -g "$MACS3_GENOMESIZE"
    -q "$MACS3_FDR"
    --format "$MACS3_FORMAT"
    --keep-dup "$MACS3_KEEPDUP"
    --outdir "$MACS3_DIR"
    --tempdir "$TMPDIR"
    --name "$name"
  )
  [[ "$use_control" -eq 1 ]] && macs_cmd+=( -c "$control_file" )
  [[ "$type_upper" == "HISTONE" ]] && macs_cmd+=( --broad --broad-cutoff "$MACS3_BROAD_CUTOFF" )

  # MACS3 does not remove a previous run's files (a row switched broad -> narrow would keep its old broadPeak, and 5_qc labels on
  # whichever file exists): start every row clean (2026-10-20)
  rm -f "${MACS3_DIR}/${name}_peaks."{narrowPeak,broadPeak,gappedPeak,xls} "${MACS3_DIR}/${name}_summits.bed" \
        "${MACS3_DIR}/${name}_filtered.bed" "${MACS3_DIR}/${name}.annotatePeaks.txt"
  "${macs_cmd[@]}"

  narrowPeak_file="${MACS3_DIR}/${name}_peaks.narrowPeak"
  broadPeak_file="${MACS3_DIR}/${name}_peaks.broadPeak"
  filtered_bed="${MACS3_DIR}/${name}_filtered.bed"

  if [[ "$SKIP_BLACKLIST" -eq 0 ]]; then
    if [[ -f "$narrowPeak_file" ]]; then
      bedtools intersect -a "$narrowPeak_file" -b "$PEAK_BLACKLIST" -v > "$filtered_bed"
    elif [[ -f "$broadPeak_file" ]]; then
      bedtools intersect -a "$broadPeak_file" -b "$PEAK_BLACKLIST" -v > "$filtered_bed"
    fi
  fi

  # Annotate peaks (filtered BED if present, else narrowPeak/broadPeak)
  bed_to_annotate=""
  if [[ -f "$filtered_bed" ]]; then
    bed_to_annotate="$filtered_bed"
  elif [[ -f "$narrowPeak_file" ]]; then
    bed_to_annotate="$narrowPeak_file"
  elif [[ -f "$broadPeak_file" ]]; then
    bed_to_annotate="$broadPeak_file"
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
