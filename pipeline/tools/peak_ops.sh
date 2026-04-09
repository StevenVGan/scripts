#!/usr/bin/env bash
# Peak set operations: intersect, distinct, or union. BED or HOMER .annotatePeaks.txt.
#
# Usage:
#   ./peak_ops.sh [OPTIONS] --mode MODE OUTPUT.bed INPUT1 INPUT2 [INPUT3 ...]
#
# Options (before positional args):
#   --mode MODE       Required. One of: intersect | distinct | union
#   --viz MODE        Post-run plots: none | upset | venn | both. Default: upset
#                     (venn uses VennDiagram for 2–4 inputs only; skipped when n>4).
#   --slop BP        Extend each peak by BP bp on both sides. Default: 0
#   --genome-sizes F Chromosome sizes file for slop (required if --slop > 0).
#   --merge          Merge overlapping regions in output (intersect/distinct only).
#   --names "N1,N2,…" Set names for UpSet / Venn (comma-separated).
#
#   OUTPUT.bed  output BED (intersect/distinct: regions in ALL; union: regions in ANY)
#   INPUT*     BED files or HOMER .annotatePeaks.txt (chr,start,end from cols 2,3,4)
#
# Example:
#   ./peak_ops.sh --mode intersect --slop 250 out.bed rep1.bed rep2.bed rep3.bed
#   ./peak_ops.sh --mode union --names "R1,R2,R3" out.bed rep1.annotatePeaks.txt rep2.bed

set -euo pipefail

SLOP_BP=0
GENOME_SIZES="/mnt/share/archive/bkup/ref/genome/hg38/hg38.chrom.sizes"
MODE=""
MERGE=0
UPSET_NAMES=""
VIZ="${VIZ:-upset}"

# Parse optional flags
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
    --merge)
      MERGE=1
      shift
      ;;
    --names)
      UPSET_NAMES="$2"
      shift 2
      ;;
    --viz)
      VIZ="$2"
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

case "$MODE" in
  intersect|distinct|union) ;;
  *)
    echo "Usage: $0 [OPTIONS] --mode intersect|distinct|union OUTPUT.bed INPUT1 INPUT2 [INPUT3 ...]" >&2
    echo "  --mode is required. Use: intersect, distinct, or union." >&2
    exit 1
    ;;
esac

case "$VIZ" in
  none|upset|venn|both) ;;
  *)
    echo "ERROR: --viz must be one of: none, upset, venn, both (got: $VIZ)" >&2
    exit 1
    ;;
esac

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 [--slop BP] [--genome-sizes FILE] [--merge] [--names \"N1,N2,...\"] --mode MODE OUTPUT.bed INPUT1 INPUT2 [INPUT3 ...]" >&2
  echo "  Requires at least 2 input files." >&2
  exit 1
fi

OUTPUT="$1"
shift
INPUTS=("$@")

command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found." >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT")"

if [[ "$SLOP_BP" -gt 0 ]]; then
  [[ -f "$GENOME_SIZES" ]] || { echo "Genome sizes file not found: $GENOME_SIZES (required for --slop)" >&2; exit 1; }
fi

# Count non-empty BED-like lines (robust vs wc -l for files without trailing newline)
count_bed_lines() {
  local c
  c=$(grep -c . "$1" 2>/dev/null || echo 0)
  c=${c%%$'\n'*}
  printf '%d' "${c:-0}"
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

# Convert all inputs to sorted BED temp files
TMPS=()
for f in "${INPUTS[@]}"; do
  T=$(mktemp)
  to_bed "$f" | sort -k1,1 -k2,2n > "$T"
  if [[ "$SLOP_BP" -gt 0 ]]; then
    T2=$(mktemp)
    bedtools slop -i "$T" -g "$GENOME_SIZES" -b "$SLOP_BP" > "$T2"
    rm -f "$T"
    T="$T2"
  fi
  TMPS+=("$T")
done

# BED generation by mode
TMP_OUT=$(mktemp)
if [[ "$MODE" == "union" ]]; then
  cat "${TMPS[@]}" | sort -k1,1 -k2,2n | bedtools merge -i stdin > "$TMP_OUT"
else
  # intersect or distinct: bedtools intersect chain
  FIRST="${TMPS[0]}"
  REST=("${TMPS[@]:1}")
  bedtools intersect -a "$FIRST" -b "${REST[0]}" > "$TMP_OUT"
  for t in "${REST[@]:1}"; do
    bedtools intersect -a "$TMP_OUT" -b "$t" > "${TMP_OUT}.n"
    mv "${TMP_OUT}.n" "$TMP_OUT"
  done
  if [[ "$MERGE" -eq 1 ]]; then
    sort -k1,1 -k2,2n "$TMP_OUT" -o "$TMP_OUT"
    bedtools merge -i "$TMP_OUT" > "${TMP_OUT}.m"
    mv "${TMP_OUT}.m" "$TMP_OUT"
  fi
fi
mv "$TMP_OUT" "$OUTPUT"

# Optional UpSet / Venn plots (--viz)
REGION_COUNT=$(count_bed_lines "$OUTPUT")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UPSET_BASE="$(dirname "$OUTPUT")/$(basename "$OUTPUT" .bed)_upset"
# peak_vennDiagram.R writes ${out_base}_venn.pdf
VENN_BASE="$(dirname "$OUTPUT")/$(basename "$OUTPUT" .bed)"
TMP_UP=$(mktemp -d)
BED_ARGS=()
for i in "${!TMPS[@]}"; do
  name=$(basename "${INPUTS[$i]}" .annotatePeaks.txt | sed 's/\.\(bed\|narrowPeak\|broadPeak\)$//')
  cp "${TMPS[$i]}" "$TMP_UP/${name}.bed"
  BED_ARGS+=("$TMP_UP/${name}.bed")
done

run_upset() {
  if [[ ! -f "${SCRIPT_DIR}/peaks_ops_upsetR.R" ]] || ! command -v Rscript >/dev/null 2>&1; then
    echo "[WARN] Rscript or peaks_ops_upsetR.R not found" >&2
    return 1
  fi
  local R_ARGS=(--mode "$MODE")
  [[ -n "$UPSET_NAMES" ]] && R_ARGS+=(--names "$UPSET_NAMES")
  Rscript "${SCRIPT_DIR}/peaks_ops_upsetR.R" "${R_ARGS[@]}" "$UPSET_BASE" "${BED_ARGS[@]}"
}

# 0 = ok (PDF written); 1 = error; 2 = skipped (n>4)
run_venn() {
  if [[ ! -f "${SCRIPT_DIR}/peak_vennDiagram.R" ]] || ! command -v Rscript >/dev/null 2>&1; then
    echo "[WARN] Rscript or peak_vennDiagram.R not found" >&2
    return 1
  fi
  if [[ "${#INPUTS[@]}" -gt 4 ]]; then
    echo "[INFO] VennDiagram supports at most 4 sets (${#INPUTS[@]} inputs); skipping Venn (--viz $VIZ)." >&2
    return 2
  fi
  local R_ARGS=()
  [[ -n "$UPSET_NAMES" ]] && R_ARGS+=(--names "$UPSET_NAMES")
  if [[ ${#R_ARGS[@]} -gt 0 ]]; then
    Rscript "${SCRIPT_DIR}/peak_vennDiagram.R" "${R_ARGS[@]}" "$VENN_BASE" "${BED_ARGS[@]}"
  else
    Rscript "${SCRIPT_DIR}/peak_vennDiagram.R" "$VENN_BASE" "${BED_ARGS[@]}"
  fi
}

case "$VIZ" in
  none) ;;
  upset)
    if run_upset; then
      echo "UpSet plot -> ${UPSET_BASE}.pdf"
    else
      echo "[WARN] peaks_ops_upsetR failed. Install: install.packages('UpSetR', repos='https://cloud.r-project.org')" >&2
    fi
    ;;
  venn)
    venn_ec=0
    run_venn || venn_ec=$?
    if [[ "$venn_ec" -eq 0 ]]; then
      echo "Venn diagram -> ${VENN_BASE}_venn.pdf"
    elif [[ "$venn_ec" -ne 2 ]]; then
      echo "[WARN] peak_vennDiagram.R failed. Install: install.packages('VennDiagram')" >&2
    fi
    ;;
  both)
    if run_upset; then
      echo "UpSet plot -> ${UPSET_BASE}.pdf"
    else
      echo "[WARN] peaks_ops_upsetR failed. Install: install.packages('UpSetR', repos='https://cloud.r-project.org')" >&2
    fi
    venn_ec=0
    run_venn || venn_ec=$?
    if [[ "$venn_ec" -eq 0 ]]; then
      echo "Venn diagram -> ${VENN_BASE}_venn.pdf"
    elif [[ "$venn_ec" -ne 2 ]]; then
      echo "[WARN] peak_vennDiagram.R failed. Install: install.packages('VennDiagram')" >&2
    fi
    ;;
esac

rm -rf "$TMP_UP"

rm -f "${TMPS[@]}"
echo "Mode=$MODE: ${#INPUTS[@]} file(s) -> $OUTPUT ($REGION_COUNT regions)"
