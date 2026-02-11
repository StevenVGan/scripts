#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 6_heatmap.sh — Heatmap of tracks vs peak/gene regions
#
# Purpose
#   Computes deepTools matrix and plots heatmaps from bigWig tracks over BED
#   regions (peaks, TSS, or gene bodies). Config-driven; no interactive prompts.
#
# Usage
#   Edit the config block below, then: ./6_heatmap.sh
#
# Output
#   OUTPUT_DIR/heatmap/, OUTPUT_DIR/sort/, OUTPUT_DIR/matrix/
#   Set KEEP_MATRIX=0 to delete matrix and sorted BED after plotting.
###############################################################################


###############################################################################
# ---- CONFIG (edit this block) ----
###############################################################################

# Base directory for this project
BASE="${HOME}/work/seq/CUTRUN/test"

# Peak sources: relative paths under PEAK_DIR, or absolute paths for external files
PEAK_DIR="${BASE}/peaks"
BED_FILES=(
  "homer/SG347_CnR_H3K4me3_MCF7_Veh_E2_vs_IgG.annotatePeaks.txt"
)
BED_LABELS=("H3K4me3")

# Track sources: relative paths under TRACK_DIR, or absolute paths
TRACK_DIR="${BASE}/align/track"
BW_FILES=(
  "SG339_CnR_ERa_MCF7_Veh_E2_rep1.bw"
  "SG347_CnR_H3K4me3_MCF7_Veh_E2.bw"
  "SG348_CnR_IgG_MCF7_Veh_ICI.bw"
)
BW_LABELS=("ERa" "H3K4me3" "IgG")

# TSS / Genes (for REGION=TSS or Genes)
TSS_BED="${HOME}/work/ref/annot/hg38/hg38_annot_genome.bed"
GENES_BED="${HOME}/work/ref/annot/hg38/hg38_annot_genome.bed"
TSS_UP=1000
TSS_DN=1000
GENES_UP=3000
GENES_DN=3000
BODY_LENGTH=7500

# Region type and peak-window params
REGION="TSS"                    # Peaks | TSS | Genes
PEAK_UP=3000
PEAK_DN=3000
OUTPUT_NAME="ERa_H3K4me3"

# Output: plot/ with subdirs heatmap/, sort/, matrix/
OUTPUT_DIR="${BASE}/plot"
LOG_DIR="${BASE}/logs"

# Plot and run options
HEATMAP_COLOR="Blues"
HEATMAP_SORT="descend"
HEATMAP_ZMIN="0"
HEATMAP_ZMAX="50"
KEEP_MATRIX=0                    # 1=keep .mtx.gz and sorted.bed; 0=delete after plot
CPU=8
# Optional: if deepTools not on PATH, set e.g. DEEPTOOLS_PREFIX="/opt/apps/bio/deeptools-3.5.1/bin"
DEEPTOOLS_PREFIX="${DEEPTOOLS_PREFIX:-}"

###############################################################################
# Helpers (local to this script)
###############################################################################

check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
}

# No need to log since we're not running this as part of the pipeline

# log_start() {
#   local step="$1"
#   mkdir -p "$LOG_DIR"
#   local logfile="${LOG_DIR}/${step}_$(date +%Y%m%d_%H%M%S).log"
#   echo "Logging to: $logfile"
#   exec > >(tee -a "$logfile") 2>&1
# }

###############################################################################
# Main
###############################################################################

# log_start "heatmap"

# Resolve deepTools commands (use DEEPTOOLS_PREFIX if set)
if [[ -n "${DEEPTOOLS_PREFIX:-}" ]]; then
  prefix="${DEEPTOOLS_PREFIX%/}/"
  COMPUTEMATRIX="${prefix}computeMatrix"
  PLOTHEATMAP="${prefix}plotHeatmap"
  [[ -x "$COMPUTEMATRIX" ]] || { echo "ERROR: $COMPUTEMATRIX not found or not executable" >&2; exit 1; }
  [[ -x "$PLOTHEATMAP" ]] || { echo "ERROR: $PLOTHEATMAP not found or not executable" >&2; exit 1; }
else
  COMPUTEMATRIX="computeMatrix"
  PLOTHEATMAP="plotHeatmap"
  check_cmd computeMatrix
  check_cmd plotHeatmap
fi

mkdir -p "${OUTPUT_DIR}/heatmap" "${OUTPUT_DIR}/sort" "${OUTPUT_DIR}/matrix"

echo "=== 6_heatmap.sh ==="
echo "[CONFIG] BASE:        $BASE"
echo "[CONFIG] REGION:      $REGION"
echo "[CONFIG] OUTPUT_NAME: $OUTPUT_NAME"
echo "[CONFIG] OUTPUT_DIR:  $OUTPUT_DIR"
echo "[CONFIG] KEEP_MATRIX: $KEEP_MATRIX"

# Resolve BED paths: relative -> under PEAK_DIR, absolute -> as-is
resolve_path() {
  local base="$1" path="$2"
  if [[ "$path" == /* ]]; then
    echo "$path"
  else
    echo "${base}/${path}"
  fi
}

# Build resolved BED list (convert HOMER annotatePeaks to BED if needed)
RESOLVED_BEDS=()
TEMP_BEDS=()
for rel in "${BED_FILES[@]}"; do
  f=$(resolve_path "$PEAK_DIR" "$rel")
  if [[ ! -f "$f" ]]; then
    echo "ERROR: BED file not found: $f" >&2
    exit 1
  fi
  if [[ "$f" == *.annotatePeaks.txt ]]; then
    # HOMER annotatePeaks: cols 2=Chr, 3=Start, 4=End
    tmpbed=$(mktemp --suffix=.bed)
    awk 'NR>1 {OFS="\t"; print $2,$3,$4}' "$f" > "$tmpbed"
    RESOLVED_BEDS+=("$tmpbed")
    TEMP_BEDS+=("$tmpbed")
  else
    RESOLVED_BEDS+=("$f")
  fi
done

# Resolve bigWig paths and check
BW_RESOLVED=()
for rel in "${BW_FILES[@]}"; do
  f=$(resolve_path "$TRACK_DIR" "$rel")
  if [[ ! -f "$f" ]]; then
    echo "ERROR: bigWig file not found: $f" >&2
    exit 1
  fi
  BW_RESOLVED+=("$f")
done
BW_FILES=("${BW_RESOLVED[@]}")

# Build space-separated strings for computeMatrix
bw_list="${BW_FILES[*]}"
bw_labels_str="${BW_LABELS[*]}"
bed_list="${RESOLVED_BEDS[*]}"
bed_labels_str="${BED_LABELS[*]}"

output_mtx="${OUTPUT_DIR}/matrix/${OUTPUT_NAME}_mtx.gz"
output_sort="${OUTPUT_DIR}/sort/${OUTPUT_NAME}_sorted.bed"
output_heatmap="${OUTPUT_DIR}/heatmap/${OUTPUT_NAME}_heatmap.png"

REGIONSLABEL_ARGS=()
if [[ "$REGION" == "Peaks" ]]; then
  REGIONSLABEL_ARGS=(--regionsLabel ${bed_labels_str})
fi

###############################################################################
# computeMatrix
###############################################################################

if [[ "$REGION" == "Peaks" ]]; then
  echo "[computeMatrix] Peaks (reference-point center) ..."
  "$COMPUTEMATRIX" reference-point --referencePoint "center" -b "$PEAK_UP" -a "$PEAK_DN" \
    -S $bw_list -R $bed_list -o "$output_mtx" --outFileSortedRegions "$output_sort" \
    --missingDataAsZero -p "$CPU" --samplesLabel $bw_labels_str

elif [[ "$REGION" == "TSS" ]]; then
  if [[ ! -f "$TSS_BED" ]]; then
    echo "ERROR: TSS_BED not found: $TSS_BED" >&2
    exit 1
  fi
  output_mtx="${OUTPUT_DIR}/matrix/${OUTPUT_NAME}_TSS_mtx.gz"
  output_sort="${OUTPUT_DIR}/sort/${OUTPUT_NAME}_TSS_sorted.bed"
  output_heatmap="${OUTPUT_DIR}/heatmap/${OUTPUT_NAME}_TSS_heatmap.png"
  echo "[computeMatrix] TSS (reference-point TSS) ..."
  "$COMPUTEMATRIX" reference-point --referencePoint "TSS" -b "$TSS_UP" -a "$TSS_DN" \
    -S $bw_list -R "$TSS_BED" -o "$output_mtx" --outFileSortedRegions "$output_sort" \
    --missingDataAsZero -p "$CPU" --samplesLabel $bw_labels_str

elif [[ "$REGION" == "Genes" ]]; then
  if [[ ! -f "$GENES_BED" ]]; then
    echo "ERROR: GENES_BED not found: $GENES_BED" >&2
    exit 1
  fi
  output_mtx="${OUTPUT_DIR}/matrix/${OUTPUT_NAME}_Genes_mtx.gz"
  output_sort="${OUTPUT_DIR}/sort/${OUTPUT_NAME}_Genes_sorted.bed"
  output_heatmap="${OUTPUT_DIR}/heatmap/${OUTPUT_NAME}_Genes_heatmap.png"
  echo "[computeMatrix] Genes (scale-regions) ..."
  "$COMPUTEMATRIX" scale-regions -m "$BODY_LENGTH" -b "$GENES_UP" -a "$GENES_DN" \
    --startLabel "TSS" --endLabel "TES" \
    -S $bw_list -R "$GENES_BED" -o "$output_mtx" --outFileSortedRegions "$output_sort" \
    --missingDataAsZero -p "$CPU" --samplesLabel $bw_labels_str

else
  echo "ERROR: REGION must be Peaks, TSS, or Genes. Got: $REGION" >&2
  exit 1
fi

echo "[computeMatrix] done -> $output_mtx"

###############################################################################
# plotHeatmap
###############################################################################

echo "[plotHeatmap] $output_mtx -> $output_heatmap"
"$PLOTHEATMAP" -m "$output_mtx" -out "$output_heatmap" \
  --colorMap "$HEATMAP_COLOR" --zMin "$HEATMAP_ZMIN" --zMax "$HEATMAP_ZMAX" \
  --missingDataColor 1 ${REGIONSLABEL_ARGS[@]+"${REGIONSLABEL_ARGS[@]}"} --sortRegions "$HEATMAP_SORT"

echo "[plotHeatmap] done"

###############################################################################
# Cleanup (if KEEP_MATRIX=0)
###############################################################################

if [[ "$KEEP_MATRIX" -eq 0 ]]; then
  echo "[cleanup] Removing matrix and sorted BED (KEEP_MATRIX=0)"
  rm -f "$output_mtx" "$output_sort"
fi

# Remove temp BED files (from HOMER annotatePeaks conversion)
for f in "${TEMP_BEDS[@]:-}"; do
  [[ -n "$f" && -f "$f" ]] && rm -f "$f"
done

echo "=== 6_heatmap.sh finished ==="
