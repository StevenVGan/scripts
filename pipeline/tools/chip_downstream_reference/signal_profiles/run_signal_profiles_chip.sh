#!/usr/bin/env bash
# ChIP-seq signal profiles / heatmaps (HOMER peaks) across 4 conditions.
# Requires: conda activate bio (bedtools, deepTools, python3).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TRACK_DIR="${TRACK_DIR:-$ROOT/../../align/track}"

OUT_DIR="${OUT_DIR:-$SCRIPT_DIR}"
GENOME_SIZES="${GENOME_SIZES:-/mnt/share/archive/bkup/ref/genome/hg38/hg38.chrom.sizes}"
PROX="${PROX:-2000}"
SAMPLE_FRAC="${SAMPLE_FRAC:-0.10}"
SEED="${SEED:-42}"
FLANK_BP="${FLANK_BP:-5000}"
NPROC="${NPROC:-max}"

NORM_MODE="${NORM_MODE:-p99}"

# Heatmap color cap and summary-profile y-axis (same units as matrix values after normalization).
# ER ChIP is much stronger than p65 here; override with ZMAX_ER, YMAX_ER, ZMAX_P65, YMAX_P65 if needed.
ZMAX_ER="${ZMAX_ER:-75}"
YMAX_ER="${YMAX_ER:-75}"
ZMAX_P65="${ZMAX_P65:-15}"
YMAX_P65="${YMAX_P65:-15}"

mkdir -p "$OUT_DIR"

for tool in computeMatrix plotHeatmap; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Missing $tool (conda activate bio?)" >&2; exit 1; }
done
command -v python3 >/dev/null 2>&1 || { echo "Missing python3" >&2; exit 1; }

declare -a CONDS=(Veh E2 TNF E2TNF)

MANIFEST="$OUT_DIR/tracks_manifest.tsv"
echo -e "factor\tcondition\trep1\trep2\tmean_bw" >"$MANIFEST"

for factor in ER p65; do
  for cond in "${CONDS[@]}"; do
    r1="${TRACK_DIR}/${factor}_${cond}_rep1.bw"
    r2="${TRACK_DIR}/${factor}_${cond}_rep2.bw"
    # We average reps at the computeMatrix stage (matrix-level mean), since bigwigCompare mean can be slow/hang on some systems.
    echo -e "${factor}\t${cond}\t${r1}\t${r2}\tMATRIX_AVG(rep1,rep2)" >>"$MANIFEST"
    [[ -f "$r1" ]] || { echo "Missing: $r1" >&2; exit 1; }
    [[ -f "$r2" ]] || { echo "Missing: $r2" >&2; exit 1; }
  done
done

echo "== Build site BEDs (HOMER-based) -> $OUT_DIR =="
python3 "$SCRIPT_DIR/build_site_beds_chip_homer.py" \
  --cobinding-root "$ROOT" \
  --genome-sizes "$GENOME_SIZES" \
  --prox "$PROX" \
  --sample-frac "$SAMPLE_FRAC" \
  --seed "$SEED" \
  --out-dir "$OUT_DIR"

COBIND="$OUT_DIR/cobind_E2TNF.bed"
P65_ONLY="$OUT_DIR/p65_TNF_only.bed"
ER_ONLY="$OUT_DIR/ER_E2_only_10pct.bed"

for f in "$COBIND" "$P65_ONLY" "$ER_ONLY"; do
  [[ -s "$f" ]] || echo "[WARN] $f is missing or empty; computeMatrix may fail." >&2
done

matrix_one_factor() {
  local factor="$1"
  local zmax="$2"
  local ymax="$3"
  local matrix_raw="$OUT_DIR/matrix_${factor}_by_siteclass.pre_norm.gz"
  local matrix="$OUT_DIR/matrix_${factor}_by_siteclass.gz"
  local fig="$OUT_DIR/figure_${factor}_by_siteclass.png"
  local bw_veh_r1="${TRACK_DIR}/${factor}_Veh_rep1.bw"
  local bw_veh_r2="${TRACK_DIR}/${factor}_Veh_rep2.bw"
  local bw_e2_r1="${TRACK_DIR}/${factor}_E2_rep1.bw"
  local bw_e2_r2="${TRACK_DIR}/${factor}_E2_rep2.bw"
  local bw_tnf_r1="${TRACK_DIR}/${factor}_TNF_rep1.bw"
  local bw_tnf_r2="${TRACK_DIR}/${factor}_TNF_rep2.bw"
  local bw_et_r1="${TRACK_DIR}/${factor}_E2TNF_rep1.bw"
  local bw_et_r2="${TRACK_DIR}/${factor}_E2TNF_rep2.bw"

  echo "== computeMatrix ($factor) -> $matrix_raw =="
  computeMatrix reference-point \
    -S "$bw_veh_r1" "$bw_veh_r2" "$bw_e2_r1" "$bw_e2_r2" "$bw_tnf_r1" "$bw_tnf_r2" "$bw_et_r1" "$bw_et_r2" \
    -R "$COBIND" "$P65_ONLY" "$ER_ONLY" \
    -b "$FLANK_BP" -a "$FLANK_BP" \
    --referencePoint center \
    --samplesLabel Veh_rep1 Veh_rep2 E2_rep1 E2_rep2 TNF_rep1 TNF_rep2 E2TNF_rep1 E2TNF_rep2 \
    -o "$matrix_raw" \
    -p "$NPROC"

  echo "== Average rep1/rep2 pairs (Veh/E2/TNF/E2TNF) -> $matrix_raw.avg.gz =="
  python3 "$SCRIPT_DIR/average_matrix_rep_pairs.py" \
    -i "$matrix_raw" \
    -o "$matrix_raw.avg.gz" \
    --pairs "Veh,Veh_rep1,Veh_rep2" "E2,E2_rep1,E2_rep2" "TNF,TNF_rep1,TNF_rep2" "E2TNF,E2TNF_rep1,E2TNF_rep2"

  if [[ "$NORM_MODE" == "none" ]]; then
    cp -f "$matrix_raw.avg.gz" "$matrix"
  else
    echo "== normalize_matrix_samples.py ($NORM_MODE) -> $matrix =="
    python3 "$SCRIPT_DIR/normalize_matrix_samples.py" \
      -i "$matrix_raw.avg.gz" \
      -o "$matrix" \
      --mode "$NORM_MODE"
  fi

  echo "== plotHeatmap ($factor, zMax=$zmax yMax=$ymax) -> $fig =="
  plotHeatmap \
    -m "$matrix" \
    -o "$fig" \
    --colorMap Blues \
    --whatToShow "plot, heatmap and colorbar" \
    --plotType lines \
    --samplesLabel Veh E2 TNF E2TNF \
    --regionsLabel Cotreat_cobind p65_only ER_only \
    --refPointLabel center \
    --yAxisLabel "Norm. signal" \
    --zMin 0 --zMax "$zmax" \
    --yMin 0 --yMax "$ymax" \
    --dpi 200
}

matrix_one_factor ER "$ZMAX_ER" "$YMAX_ER"
matrix_one_factor p65 "$ZMAX_P65" "$YMAX_P65"

echo "Done. Outputs under: $OUT_DIR"
