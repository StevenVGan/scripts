#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# spikein_ecoli.sh — E. coli spike-in readout for a cnr project (CUT&RUN / CUT&Tag kits that add E. coli DNA)
#
# Promoted 2026-11-01 from seq/cnr/CnR_261006_ERa-ARNT-p65_MCF7-HEK293T_Steven/script/local/spikein_ecoli.sh
# (CONVENTIONS §5: needed in a third project). Behaviour unchanged; SPIKEIN_QC_DIR and SPIKEIN_MULTIQC added so a
# side-by-side re-run (script_r2/, cleandata_r2/) can write its own table without touching the project's qc/ or report.
#
# Purpose
#   Aligns the trimmed pairs to E. coli K-12 MG1655 and tabulates, per sample, hg38 vs E. coli concordant pairs.
#   With a fixed spike-in, hg38 / E. coli is the recovery readout (Qubit yields cannot judge TF arms).
#
# Inputs (from the project config: CONFIG_FILE, default ./0_config.sh, i.e. run from the project's script/ dir)
#   - Trimmed pairs ${TRIM_DIR}/*_R1_val_1.fq.gz + *_R2_val_2.fq.gz (set DELETE_TRIMMED_AFTER_ALIGN=0 so step 2 keeps them)
#   - hg38 BAMs ${BAM_DIR}/<sample>_sorted.bam (step 2)
#   - E. coli bowtie2 index ${ECOLI_INDEX} (NC_000913.3 K-12 MG1655, 4,641,652 bp)
#
# Outputs (SPIKEIN_QC_DIR, default ${BASE}/qc)
#   - spikein/<sample>_ecoli_bowtie2.log
#   - spikein_ecoli.tsv  one row per sample:
#       sample  input_pairs  hg38_pairs  ecoli_pairs  pct_ecoli  hg38_per_ecoli  scale_factor
#     scale_factor = SCALE_CONST / ecoli_pairs (Henikoff-style; not applied to any track here)
#   - spikein_ecoli_gs_mqc.tsv + spikein_ecoli_table_mqc.tsv (MultiQC custom content: 2 General Statistics columns +
#     an "E. coli spike-in" table); with SPIKEIN_MULTIQC=1 (default) MultiQC is re-run with 5_qc.sh's config so the
#     report carries them (*_ecoli_bowtie2.log ignored — they would otherwise land in the Bowtie2 section).
#
# Behavior
#   - Re-aligns a sample only if its E. coli log is missing or older than its trimmed R1.
#   - Pairing as the cnr default hg38 call (--dovetail, concordant only, -I 10 -X 700). The Henikoff recipe's
#     --no-overlap was written for PE25: at PE101 it makes every fragment < 202 bp non-concordant.
#   - DELETE_TRIMMED_AFTER_SPIKEIN=1 deletes the trimmed FASTQs once EVERY sample has a row.
#
# Usage, after run_all.sh (needs step 2's BAMs and the kept trimmed FASTQs; step 5's multiqc_config.yaml for the report):
#   cd <project>/script && ~/work/scripts/pipeline/tools/spikein_ecoli.sh
#   CONFIG_FILE=<project>/script_r2/0_config.sh SPIKEIN_QC_DIR=<project>/qc_r2 SPIKEIN_MULTIQC=0 ~/work/scripts/pipeline/tools/spikein_ecoli.sh
###############################################################################

CONFIG_FILE="${CONFIG_FILE:-${PWD}/0_config.sh}"
[[ -f "$CONFIG_FILE" ]] || { echo "[ERROR] project config not found: $CONFIG_FILE (run from the project's script/ dir or set CONFIG_FILE)" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG_FILE"

# ---- config ----
ECOLI_INDEX="${ECOLI_INDEX:-$HOME/work/ref/align/bowtie2/cellculture_v1/ecoli_NC_000913.3/ecoli_NC_000913.3}"   # NC_000913.3, same record as the contaminant panel (2026-10-20)
SPIKE_CPU="${SPIKE_CPU:-12}"
SCALE_CONST="${SCALE_CONST:-10000}"
DELETE_TRIMMED_AFTER_SPIKEIN="${DELETE_TRIMMED_AFTER_SPIKEIN:-0}"   # 1 = delete the trimmed FASTQs once every sample has a row
SPIKEIN_QC_DIR="${SPIKEIN_QC_DIR:-${BASE}/qc}"
SPIKEIN_MULTIQC="${SPIKEIN_MULTIQC:-1}"
SPIKE_DIR="${SPIKEIN_QC_DIR}/spikein"
OUT_TSV="${SPIKEIN_QC_DIR}/spikein_ecoli.tsv"

log_start "spikein_ecoli"
check_cmd bowtie2
check_cmd samtools
mkdir -p "$SPIKE_DIR"
compgen -G "${ECOLI_INDEX}.1.bt2*" > /dev/null || { echo "[ERROR] E. coli bowtie2 index not found: ${ECOLI_INDEX}" >&2; exit 1; }

echo "=== E. coli spike-in ==="
echo "[INFO] CONFIG_FILE:    $CONFIG_FILE"
echo "[INFO] TRIM_DIR:       $TRIM_DIR"
echo "[INFO] BAM_DIR:        $BAM_DIR"
echo "[INFO] SPIKEIN_QC_DIR: $SPIKEIN_QC_DIR"
echo "[INFO] ECOLI_INDEX:    $ECOLI_INDEX"

r1_files=( "${TRIM_DIR}"/*_R1_val_1.fq.gz )
(( ${#r1_files[@]} > 0 )) || { echo "[ERROR] no *_R1_val_1.fq.gz in $TRIM_DIR (deleted? rerun step 1 with DELETE_TRIMMED_AFTER_ALIGN=0)" >&2; exit 1; }

printf 'sample\tinput_pairs\thg38_pairs\tecoli_pairs\tpct_ecoli\thg38_per_ecoli\tscale_factor\n' > "${OUT_TSV}.tmp"
n_ok=0
for R1 in "${r1_files[@]}"; do
  R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"
  sample="$(basename "$R1" _R1_val_1.fq.gz)"
  [[ -f "$R2" ]] || { echo "[ERROR] missing R2 for $sample" >&2; exit 1; }
  bam="${BAM_DIR}/${sample}_sorted.bam"
  [[ -f "$bam" ]] || { echo "[ERROR] hg38 BAM missing for $sample: $bam (run step 2 first)" >&2; exit 1; }
  eco_log="${SPIKE_DIR}/${sample}_ecoli_bowtie2.log"

  if [[ -s "$eco_log" && "$eco_log" -nt "$R1" ]]; then
    echo "[SPIKE] $sample: E. coli log up to date, reusing"
  else
    echo "[SPIKE] $sample: aligning to E. coli"
    bowtie2 -x "$ECOLI_INDEX" -1 "$R1" -2 "$R2" --end-to-end --very-sensitive --dovetail \
      --no-mixed --no-discordant --no-unal -I 10 -X 700 -p "$SPIKE_CPU" 2> "${eco_log}.tmp" > /dev/null
    mv -f "${eco_log}.tmp" "$eco_log"
  fi

  input_pairs=$(awk '/reads; of these:/{print $1; exit}' "$eco_log")
  # concordant pairs = aligned concordantly exactly 1 time + >1 times
  ecoli_pairs=$(awk '/aligned concordantly exactly 1 time|aligned concordantly >1 times/{s+=$1} END{print s+0}' "$eco_log")
  hg38_pairs=$(( $(samtools view -c -f 2 -F 0x904 -@ 4 "$bam") / 2 ))
  [[ -n "$input_pairs" ]] || { echo "[ERROR] could not parse $eco_log" >&2; exit 1; }

  awk -v s="$sample" -v i="$input_pairs" -v h="$hg38_pairs" -v e="$ecoli_pairs" -v c="$SCALE_CONST" 'BEGIN{
    printf "%s\t%d\t%d\t%d\t%.3f\t%s\t%s\n", s, i, h, e, (i>0?100*e/i:0),
      (e>0?sprintf("%.2f",h/e):"NA"), (e>0?sprintf("%.6g",c/e):"NA") }' >> "${OUT_TSV}.tmp"
  n_ok=$((n_ok + 1))
done
mv -f "${OUT_TSV}.tmp" "$OUT_TSV"
echo "[SPIKE] wrote $OUT_TSV ($n_ok samples)"
column -t -s $'\t' "$OUT_TSV"

# ---- MultiQC custom content ----
{
  cat <<'HDR'
# id: 'spikein_ecoli_gs'
# plot_type: 'generalstats'
# pconfig:
#     pct_ecoli:
#         title: '% E. coli'
#         description: 'Concordant pairs aligned to E. coli K-12 MG1655, % of trimmed input pairs'
#         min: 0
#         suffix: '%'
#         format: '{:,.2f}'
#     hg38_per_ecoli:
#         title: 'hg38/E. coli'
#         description: 'hg38 properly paired pairs per E. coli pair (fixed spike-in: higher = more target DNA)'
#         format: '{:,.1f}'
HDR
  printf 'Sample\tpct_ecoli\thg38_per_ecoli\n'
  awk -F'\t' 'NR>1{print $1"\t"$5"\t"$6}' "$OUT_TSV"
} > "${SPIKEIN_QC_DIR}/spikein_ecoli_gs_mqc.tsv"
{
  cat <<'HDR'
# id: 'spikein_ecoli_table'
# section_name: 'E. coli spike-in'
# description: 'hg38 vs E. coli concordant pairs per sample (scripts/pipeline/tools/spikein_ecoli.sh). A defined E. coli spike-in went into every reaction, so hg38/E. coli, not Qubit yield, is the recovery readout. scale_factor = 10000 / ecoli_pairs (not applied to tracks).'
# plot_type: 'table'
HDR
  sed '1s/^sample/Sample/' "$OUT_TSV"
} > "${SPIKEIN_QC_DIR}/spikein_ecoli_table_mqc.tsv"
if [[ "$SPIKEIN_MULTIQC" -ne 1 ]]; then
  echo "[SPIKE] SPIKEIN_MULTIQC=0: report not refreshed"
elif [[ -f "${BASE}/multiqc_config.yaml" ]]; then
  check_cmd multiqc
  multiqc "$BASE" "$BAMQC_DIR" -o "$MULTIQC_OUT" -c "${BASE}/multiqc_config.yaml" -f --ignore "*_ecoli_bowtie2.log"
  echo "[SPIKE] MultiQC report refreshed with the spike-in columns"
else
  echo "[WARN] ${BASE}/multiqc_config.yaml missing (5_qc.sh not run?) — spike-in _mqc files written, report NOT refreshed" >&2
fi

if [[ "$DELETE_TRIMMED_AFTER_SPIKEIN" -eq 1 && "$n_ok" -eq "${#r1_files[@]}" ]]; then
  rm -f "${TRIM_DIR}"/*_val_*.fq.gz
  echo "[SPIKE] trimmed FASTQs deleted (DELETE_TRIMMED_AFTER_SPIKEIN=1)"
fi
echo "=== E. coli spike-in complete ==="
