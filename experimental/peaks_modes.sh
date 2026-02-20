#!/usr/bin/env bash
# Peak set operations with ComplexHeatmap UpSet (distinct / intersect / union modes).
#
# Purpose
#   Single command to run peak set analysis in one of three modes. Each mode produces
#   an UpSet PDF and a BED file of regions for the "all" combination.
#
# Usage:
#   ./peaks_modes.sh [OPTIONS] --mode MODE OUTPUT_BASE INPUT1 INPUT2 [INPUT3 ...]
#
# Options (before positional args):
#   --mode MODE       Required. One of: distinct | intersect | union
#     distinct:  Mutually exclusive partitions. BED = regions in ALL inputs.
#     intersect: 1=in, 0=ignored. BED = regions in ALL inputs.
#     union:     1=in, 0=ignored, OR logic. BED = regions in ANY input.
#   --slop BP         Extend peaks by BP bp before analysis. Default: 0
#   --genome-sizes F  Chromosome sizes for slop. Default: hg38.chrom.sizes
#   --names "N1,N2,…" Set names for UpSet (comma-separated).
#
# Output:
#   OUTPUT_BASE.pdf   UpSet plot
#   OUTPUT_BASE.bed   Regions: intersect of all (distinct/intersect) or union of all (union)
#
# Example:
#   ./peaks_modes.sh --mode intersect --slop 250 --names "R1,R2,R3,R4" peaks/ERa_all peaks/*.bed
#
# Requires: R + ComplexHeatmap + GenomicRanges
#   Install: R -e "BiocManager::install(c('ComplexHeatmap','GenomicRanges'))"

set -euo pipefail

MODE=""
SLOP_BP=0
GENOME_SIZES="/mnt/share/archive/bkup/ref/genome/hg38/hg38.chrom.sizes"
NAMES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    --slop)
      SLOP_BP="$2"
      shift 2
      ;;
    --genome-sizes)
      GENOME_SIZES="$2"
      shift 2
      ;;
    --names)
      NAMES="$2"
      shift 2
      ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

if [[ -z "$MODE" ]]; then
  echo "Usage: $0 [--slop BP] [--genome-sizes F] [--names \"N1,N2,...\"] --mode distinct|intersect|union OUTPUT_BASE INPUT1 INPUT2 [INPUT3 ...]" >&2
  echo "  --mode is required." >&2
  exit 1
fi

case "$MODE" in
  distinct|intersect|union) ;;
  *)
    echo "Invalid --mode: $MODE. Use distinct, intersect, or union." >&2
    exit 1
    ;;
esac

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 [OPTIONS] --mode MODE OUTPUT_BASE INPUT1 INPUT2 [INPUT3 ...]" >&2
  echo "  Requires at least 2 input files." >&2
  exit 1
fi

OUTPUT_BASE="$1"
shift
INPUTS=("$@")

command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found." >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found." >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT_BASE")"

[[ "$SLOP_BP" -gt 0 ]] && [[ ! -f "$GENOME_SIZES" ]] && {
  echo "Genome sizes not found: $GENOME_SIZES (required for --slop)" >&2
  exit 1
}

to_bed() {
  local f="$1"
  if [[ "$f" == *.annotatePeaks.txt ]]; then
    awk 'NR>1 {OFS="\t"; print $2,$3,$4}' "$f"
  else
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$f"
  fi
}

for f in "${INPUTS[@]}"; do
  [[ -f "$f" ]] || { echo "Not found: $f" >&2; exit 1; }
done

# Convert inputs to BED (with slop)
TMP_UP=$(mktemp -d)
BED_ARGS=()
for i in "${!INPUTS[@]}"; do
  name=$(basename "${INPUTS[$i]}" .annotatePeaks.txt)
  out="$TMP_UP/${name}.bed"
  to_bed "${INPUTS[$i]}" | sort -k1,1 -k2,2n > "$out"
  if [[ "$SLOP_BP" -gt 0 ]]; then
    bedtools slop -i "$out" -g "$GENOME_SIZES" -b "$SLOP_BP" > "${out}.s"
    mv "${out}.s" "$out"
  fi
  BED_ARGS+=("$out")
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_ARGS=()
[[ -n "$NAMES" ]] && R_ARGS+=(--names "$NAMES")
if Rscript "${SCRIPT_DIR}/peaks_modes_complexheatmap.R" "${R_ARGS[@]}" "$MODE" "$OUTPUT_BASE" "${BED_ARGS[@]}"; then
  echo "Done. Mode=$MODE -> ${OUTPUT_BASE}.pdf, ${OUTPUT_BASE}.bed"
else
  rm -rf "$TMP_UP"
  echo "[ERROR] R script failed. Install: BiocManager::install(c('ComplexHeatmap','GenomicRanges'))" >&2
  exit 1
fi

rm -rf "$TMP_UP"
