#!/usr/bin/env bash
# Extract sequences from reference genome for BED regions (bedtools getfasta).
#
# Purpose
#   Retrieves FASTA sequences for genomic regions. By default, centers each region
#   and extends ±slop bp (e.g. 500 bp → 1001 bp total) for motif/sequence analysis.
#
# Usage:
#   ./getfasta.sh [OPTIONS] OUTPUT.txt INPUT1.bed [INPUT2.bed ...]
#
# Options (before positional args):
#   --slop BP         Extend each region by BP bp on both sides (center-based).
#                     Default: 500 (→ 1001 bp total when centering)
#   --genome-sizes F  Chromosome sizes file for slop (required if --slop > 0).
#                     Default: /mnt/share/archive/ref/genome/hg38/hg38.chrom.sizes
#   --genome-fasta F  Reference genome FASTA.
#                     Default: /mnt/share/archive/ref/genome/hg38/hg38.fa
#   --no-center       Use original BED coordinates as-is (no centering).
#   --tab             Output tab format (name, seq) instead of sequences only.
#
#   OUTPUT.txt  output file (sequences only, or tab format if --tab)
#   INPUT*.bed  BED files (chr, start, end)
#
# Example:
#   ./getfasta.sh seqs.txt peaks.bed
#   ./getfasta.sh --slop 250 --no-center seqs.txt rep1.bed rep2.bed

set -euo pipefail

SLOP_BP=500
GENOME_SIZES="/mnt/share/archive/bkup/ref/genome/hg38/hg38.chrom.sizes"
GENOME_FASTA="/mnt/share/archive/bkup/ref/genome/hg38/hg38.fa"
CENTER=1
TAB_OUT=0

# Parse optional flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    --slop)
      SLOP_BP="$2"
      shift 2
      ;;
    --genome-sizes)
      GENOME_SIZES="$2"
      shift 2
      ;;
    --genome-fasta)
      GENOME_FASTA="$2"
      shift 2
      ;;
    --no-center)
      CENTER=0
      shift
      ;;
    --tab)
      TAB_OUT=1
      shift
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

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 [--slop BP] [--genome-sizes FILE] [--genome-fasta FILE] [--no-center] [--tab] OUTPUT.txt INPUT1.bed [INPUT2.bed ...]" >&2
  echo "  Requires at least 1 input BED file." >&2
  exit 1
fi

OUTPUT="$1"
shift
INPUTS=("$@")

command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found." >&2; exit 1; }
[[ -f "$GENOME_FASTA" ]] || { echo "Genome FASTA not found: $GENOME_FASTA" >&2; exit 1; }
if [[ "$SLOP_BP" -gt 0 ]]; then
  [[ -f "$GENOME_SIZES" ]] || { echo "Genome sizes file not found: $GENOME_SIZES (required for --slop)" >&2; exit 1; }
fi
mkdir -p "$(dirname "$OUTPUT")"

for f in "${INPUTS[@]}"; do
  [[ -f "$f" ]] || { echo "Not found: $f" >&2; exit 1; }
done

# Build BED: center regions if requested, then optionally slop
to_bed() {
  local f="$1"
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$f"
}

center_bed() {
  awk 'BEGIN{OFS="\t"} {
    center=int(($2+$3)/2)
    print $1, center, center+1
  }'
}

TMP_BED=$(mktemp)
for f in "${INPUTS[@]}"; do
  to_bed "$f"
done | sort -k1,1 -k2,2n | uniq > "$TMP_BED"

if [[ "$CENTER" -eq 1 ]]; then
  TMP_CENTER=$(mktemp)
  center_bed < "$TMP_BED" > "$TMP_CENTER"
  mv "$TMP_CENTER" "$TMP_BED"
fi

if [[ "$SLOP_BP" -gt 0 ]]; then
  TMP_SLOP=$(mktemp)
  bedtools slop -i "$TMP_BED" -g "$GENOME_SIZES" -b "$SLOP_BP" > "$TMP_SLOP"
  mv "$TMP_SLOP" "$TMP_BED"
fi

TMP_FASTA=$(mktemp)
bedtools getfasta -tab -fi "$GENOME_FASTA" -bed "$TMP_BED" -fo "$TMP_FASTA"

if [[ "$TAB_OUT" -eq 1 ]]; then
  # Uppercase sequences, keep tab format
  awk 'BEGIN{OFS="\t"} {print $1, toupper($2)}' "$TMP_FASTA" > "$OUTPUT"
else
  # Sequences only, uppercase
  awk '{print toupper($2)}' "$TMP_FASTA" > "$OUTPUT"
fi

rm -f "$TMP_BED" "$TMP_FASTA"

REGION_COUNT=$(wc -l < "$OUTPUT")
echo "getfasta: ${#INPUTS[@]} file(s) -> $OUTPUT ($REGION_COUNT sequences)"
