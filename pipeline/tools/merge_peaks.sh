#!/usr/bin/env bash
# Merge peak files (BED or HOMER .annotatePeaks.txt) into one BED.
#
# Usage:
#   ./merge_peaks.sh OUTPUT.bed INPUT1 [INPUT2 ...]
#
#   OUTPUT.bed  output merged BED
#   INPUT*      BED files or HOMER .annotatePeaks.txt (chr,start,end from cols 2,3,4)
#
# Example:
#   ./merge_peaks.sh merged.bed rep1_peaks.bed rep2.annotatePeaks.txt

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 OUTPUT.bed INPUT1 [INPUT2 ...]" >&2
  exit 1
fi

OUTPUT="$1"
shift
INPUTS=("$@")

command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found." >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT")"

to_bed() {
  local f="$1"
  if [[ "$f" == *.annotatePeaks.txt ]]; then
    awk 'NR>1 {OFS="\t"; print $2,$3,$4}' "$f"
  else
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$f"
  fi
}

TMP=$(mktemp)
for f in "${INPUTS[@]}"; do
  [[ -f "$f" ]] || { echo "Not found: $f" >&2; exit 1; }
  to_bed "$f"
done | sort -k1,1 -k2,2n > "$TMP"

bedtools merge -i "$TMP" > "$OUTPUT"
rm -f "$TMP"
echo "Merged ${#INPUTS[@]} file(s) -> $OUTPUT ($(wc -l < "$OUTPUT") regions)"
