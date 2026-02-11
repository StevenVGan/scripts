#!/usr/bin/env bash
# General BED lifter: user-defined input (dir or files), assembly FROM → TO, output dir.
#
# Usage:
#   ./lift_bed.sh INPUT FROM TO OUTPUT [CHAIN]
#
#   INPUT   directory (all *.bed inside) or a single .bed file
#   FROM    source assembly (e.g. hg19)
#   TO      target assembly (e.g. hg38)
#   OUTPUT  output directory
#   CHAIN   optional; default: script dir / ${FROM}To${TO}.over.chain
#
# Examples:
#   ./lift_bed.sh ./MCF7_Amir_hg19 hg19 hg38 ./MCF7_Amir_hg38
#   ./lift_bed.sh my.bed hg19 hg38 ./out
#
# Chain (once): wget then gunzip from
#   http://hgdownload.soe.ucsc.edu/goldenPath/<FROM>/liftOver/<FROM>To<TO>.over.chain.gz

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 INPUT FROM TO OUTPUT [CHAIN]" >&2
  echo "  INPUT:  dir (all *.bed) or single .bed file" >&2
  echo "  FROM/TO: e.g. hg19, hg38" >&2
  exit 1
fi

INPUT="$1"
FROM="$2"
TO="$3"
OUT_DIR="$4"
CHAIN="${5:-${SCRIPT_DIR}/${FROM}To${TO}.over.chain}"



if [[ ! -f "$CHAIN" ]]; then
  echo "Downloading chain file from UCSC..."
  TO_CAPITALIZED=$(echo "$TO" | sed 's/^./\U&/')
  wget -O "${CHAIN}.gz" "http://hgdownload.soe.ucsc.edu/goldenPath/${FROM}/liftOver/${FROM}To${TO_CAPITALIZED}.over.chain.gz"
  gunzip -k "${CHAIN}.gz"
else
  echo "Using existing chain file: $CHAIN"
fi

mkdir -p "$OUT_DIR"

BEDS=()
if [[ -d "$INPUT" ]]; then
  for f in "$INPUT"/*.bed; do
    [[ -f "$f" ]] && BEDS+=("$f")
  done
else
  [[ -f "$INPUT" ]] && BEDS+=("$INPUT")
fi

[[ ${#BEDS[@]} -eq 0 ]] && echo "No BED files found." >&2 && exit 1

# liftOver expects standard BED; extended BEDs (e.g. encode cCRE with text in col10) break it. Use first 4 cols only.
UNMAPPED_TMP=$(mktemp)
for bed in "${BEDS[@]}"; do
  base="$(basename "$bed" .bed)"
  liftOver <(cut -f1-4 "$bed") "$CHAIN" "${OUT_DIR}/${base}.bed" "$UNMAPPED_TMP"
  n_unmapped=$(wc -l < "$UNMAPPED_TMP")
  echo "Lifted: $base  ($n_unmapped unmapped)"
done
rm -f "$UNMAPPED_TMP"

echo "Done → $OUT_DIR"
