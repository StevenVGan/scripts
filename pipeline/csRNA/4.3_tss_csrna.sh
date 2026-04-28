#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.3_tss_csrna.sh — HOMER findcsRNATSS.pl (csRNA-specific TSS calling)
#
# Why this step exists
#   csRNA-seq captures 5′-capped short RNA ends, so "peaks" are single-nt TSS
#   stacks rather than ChIP-style signal domains. Calling these with MACS3 /
#   HOMER findPeaks mis-models both the shape and the noise (stable short
#   RNAs / degradation products are the dominant false positives). findcsRNATSS
#   compares cap-enriched reads against a short-RNA input tag directory and
#   emits annotated TSS clusters (Duttke et al. 2019, Genome Research).
#
# Input file — TSS_GROUPS_FILE (default: ${BASE}/tss_groups.tsv):
#     ip_bam  control_bam  name
#   Strict 1:1 mapping. No pooling / merging — if you need a different input
#   for a sample, edit that row. Rows with control "-"/none/NA are skipped.
#
# Outputs (under ${TSS_DIR}/<name>/):
#   <name>.tss.txt           HOMER peak file (valid TSS clusters; this is the main file)
#   <name>.alltss.txt        all candidate TSSs before input/rna filtering
#   <name>.stats.txt         run summary (counts by class, thresholds used)
#   <name>.inputDistribution.txt, <name>.tss.txt.freq.tsv  diagnostic tables
#
# To convert to BED: `pos2bed.pl <name>.tss.txt > <name>.tss.bed` (HOMER utility).
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.3_tss_csrna"

check_cmd findcsRNATSS.pl

echo "=== STEP 4.3: findcsRNATSS.pl (csRNA TSS calling) ==="
echo "[INFO] Groups file:          $TSS_GROUPS_FILE"
echo "[INFO] TAG_DIR:              $TAG_DIR"
echo "[INFO] TSS_DIR:              $TSS_DIR"
echo "[INFO] CSRNA_NTAG_THRESHOLD: $CSRNA_NTAG_THRESHOLD"
echo "[INFO] CSRNA_GTF:            ${CSRNA_GTF:-(none — using -genome $GENOME)}"

if [[ ! -f "$TSS_GROUPS_FILE" ]]; then
  echo "[ERROR] Groups file not found: $TSS_GROUPS_FILE"
  exit 1
fi

mkdir -p "$TSS_DIR"

gtf_opts=()
if [[ -n "${CSRNA_GTF:-}" ]]; then
  if [[ ! -f "$CSRNA_GTF" ]]; then
    echo "[ERROR] CSRNA_GTF set but file not found: $CSRNA_GTF"
    exit 1
  fi
  gtf_opts=( -gtf "$CSRNA_GTF" )
fi

while IFS=$'\t' read -r ip_spec ctrl_spec name _rest; do
  [[ -z "${ip_spec:-}" ]] && continue
  [[ "$ip_spec" =~ ^# ]] && continue

  if [[ -z "${name:-}" ]]; then name="$(basename "$ip_spec" .bam)"; fi

  # Skip rows without a valid control — findcsRNATSS needs -i
  if [[ -z "${ctrl_spec:-}" ]] || [[ "$ctrl_spec" == "-" ]] || [[ "$ctrl_spec" =~ ^(none|NONE|NA)$ ]]; then
    echo "[WARN] $name: no control in TSV — skipping (findcsRNATSS needs -i)"
    continue
  fi

  # Resolve IP BAM → tagdir
  if [[ -f "$ip_spec" ]]; then
    ip_file="$ip_spec"
  elif [[ -f "${BAM_DIR}/${ip_spec}" ]]; then
    ip_file="${BAM_DIR}/${ip_spec}"
  else
    echo "[WARN] IP BAM not found for '$ip_spec' — skipping $name"
    continue
  fi
  ip_tag_dir="${TAG_DIR}/$(basename "$ip_file" _sorted.bam)"
  if [[ ! -d "$ip_tag_dir" ]]; then
    echo "[WARN] IP tag dir not found: $ip_tag_dir — skipping $name"
    continue
  fi

  # Resolve control BAM → tagdir
  if [[ -f "$ctrl_spec" ]]; then
    ctrl_bam="$ctrl_spec"
  elif [[ -f "${BAM_DIR}/${ctrl_spec}" ]]; then
    ctrl_bam="${BAM_DIR}/${ctrl_spec}"
  else
    echo "[WARN] Control BAM not found for '$ctrl_spec' — skipping $name"
    continue
  fi
  ctrl_tag_dir="${TAG_DIR}/$(basename "$ctrl_bam" _sorted.bam)"
  if [[ ! -d "$ctrl_tag_dir" ]]; then
    echo "[WARN] Control tag dir not found: $ctrl_tag_dir — skipping $name"
    continue
  fi

  outdir="${TSS_DIR}/${name}"
  out_prefix="${outdir}/${name}"
  # Require both main and stats outputs (stats.txt is written at the end — its
  # presence together with a non-empty tss.txt indicates a completed run and
  # avoids treating a half-written output from a killed previous run as "done").
  if [[ -s "${out_prefix}.tss.txt" && -f "${out_prefix}.stats.txt" ]]; then
    echo "[STEP4.3] $name: skipping — ${out_prefix}.tss.txt + .stats.txt exist"
    continue
  fi

  echo "[STEP4.3] $name: IP=$(basename "$ip_tag_dir")  input=$(basename "$ctrl_tag_dir")"
  mkdir -p "$outdir"

  findcsRNATSS.pl "$ip_tag_dir" \
    -o "$out_prefix" \
    -i "$ctrl_tag_dir" \
    -genome "$GENOME" \
    -ntagThreshold "$CSRNA_NTAG_THRESHOLD" \
    ${gtf_opts[@]+"${gtf_opts[@]}"}
done < "$TSS_GROUPS_FILE"

echo "=== STEP 4.3 complete ==="
