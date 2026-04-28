#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 2_bowtie2.sh — ATAC alignment + dedup + filter + Tn5 shift + bigWig
#
# Per-sample flow (PE)
#   bowtie2 -X $BT2_MAX_FRAG --no-mixed --no-discordant --very-sensitive
#     |  samtools sort -n      (name sort for fixmate)
#     |  samtools fixmate -m   (add MS/ms tags for markdup)
#     |  samtools sort         (coord sort)
#     |  samtools markdup      (mark dups, do NOT remove yet)
#     -> ${sample}_marked.bam              (kept for preseq complexity QC)
#
#   samtools view -b -q $MAPQ_MIN -F $EXCLUDE_FLAGS -f 2 \
#     ${sample}_marked.bam <non-chrM contigs from idxstats>
#     -> ${sample}_sorted.bam               (filtered + dedup'd, autosomes/sex only)
#
#   alignmentSieve --ATACshift -b ${sample}_sorted.bam
#     -> sort -> ${sample}_sorted.shifted.bam   (Tn5 +4/-5 cut-site shift, used downstream)
#
#   bamCoverage on ${sample}_sorted.shifted.bam -> ${sample}.bw
#
# SE: drops -f 2 from the post-filter view; markdup still works (position+orientation).
#
# Idempotency: skips alignment / shift / bigWig if outputs exist.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "2_bowtie2"

check_cmd bowtie2
check_cmd samtools
check_cmd alignmentSieve
check_cmd bamCoverage

echo "=== STEP 2: bowtie2 -> dedup -> filter -> Tn5 shift -> bigWig ==="
echo "[INFO] GENOME_INDEX: $GENOME_INDEX"
echo "[INFO] BT2_MAX_FRAG: $BT2_MAX_FRAG"
echo "[INFO] MAPQ_MIN:     $MAPQ_MIN"
echo "[INFO] EXCLUDE_CHROMS: $EXCLUDE_CHROMS"
echo "[INFO] Mode:         $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"

# build a sed/grep pattern for chrM exclusion: "chrM\|chrEBV\|^\*$"
chrom_exclude_pat='^\*$'
for c in $EXCLUDE_CHROMS; do
  chrom_exclude_pat="${chrom_exclude_pat}|^${c}$"
done

# shellcheck disable=SC2206
bt2_extra_args=( $BT2_EXTRA )

# ---- helpers ----
align_and_dedup() {
  # $1=R1, $2=R2 (empty for SE), $3=marked_bam (output)
  local R1="$1" R2="$2" marked="$3" log="$4"
  if [[ "${SE:-0}" -eq 1 ]]; then
    bowtie2 -x "$GENOME_INDEX" -U "$R1" \
      --very-sensitive -X "$BT2_MAX_FRAG" -p "$BT2_CPU" \
      "${bt2_extra_args[@]}" 2> "$log" \
      | samtools sort -@ "$BT2_CPU" -o - - \
      | samtools markdup -@ "$BT2_CPU" - "$marked"
  else
    bowtie2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" \
      --very-sensitive -X "$BT2_MAX_FRAG" -p "$BT2_CPU" \
      "${bt2_extra_args[@]}" 2> "$log" \
      | samtools sort -n -@ "$BT2_CPU" - \
      | samtools fixmate -m -@ "$BT2_CPU" - - \
      | samtools sort -@ "$BT2_CPU" - \
      | samtools markdup -@ "$BT2_CPU" - "$marked"
  fi
  samtools index "$marked"
}

filter_bam() {
  # ${sample}_marked.bam -> ${sample}_sorted.bam (filtered + dedup'd, chrM removed)
  local marked="$1" sorted="$2"
  local view_args=( -b -q "$MAPQ_MIN" -F "$EXCLUDE_FLAGS" )
  [[ "${SE:-0}" -eq 0 ]] && view_args+=( -f 2 )

  # Build region list of contigs we want to keep (everything in idxstats minus excluded)
  mapfile -t keep_chroms < <(samtools idxstats "$marked" \
    | cut -f1 \
    | grep -E -v "$chrom_exclude_pat")

  if (( ${#keep_chroms[@]} == 0 )); then
    echo "[ERROR] No contigs left after chrM/* filter — check EXCLUDE_CHROMS." >&2
    exit 1
  fi

  samtools view "${view_args[@]}" -@ "$BT2_CPU" "$marked" "${keep_chroms[@]}" -o "$sorted"
  samtools index "$sorted"
}

tn5_shift() {
  # ${sample}_sorted.bam -> ${sample}_sorted.shifted.bam (sort+index)
  local sorted="$1" shifted="$2" tmp
  tmp="${shifted%.bam}.unsorted.bam"
  alignmentSieve --ATACshift --bam "$sorted" --outFile "$tmp" -p "$BT2_CPU"
  samtools sort -@ "$BT2_CPU" -o "$shifted" "$tmp"
  samtools index "$shifted"
  rm -f "$tmp"
}

# ---- main loop ----
if [[ "${SE:-0}" -eq 1 ]]; then
  trimmed_files=( "${TRIM_DIR}"/*_trimmed.fq.gz )
  if (( ${#trimmed_files[@]} == 0 )); then
    echo "[ERROR] No *_trimmed.fq.gz files found in $TRIM_DIR"
    exit 1
  fi
  inputs=( "${trimmed_files[@]}" )
else
  r1_files=( "${TRIM_DIR}"/*_R1_val_1.fq.gz )
  if (( ${#r1_files[@]} == 0 )); then
    echo "[ERROR] No *_R1_val_1.fq.gz files found in $TRIM_DIR"
    exit 1
  fi
  inputs=( "${r1_files[@]}" )
fi

for R1 in "${inputs[@]}"; do
  if [[ "${SE:-0}" -eq 1 ]]; then
    base="$(basename "$R1" .fq.gz)"
    [[ "$base" == *_R1_trimmed ]] && sample_name="${base%_R1_trimmed}" || sample_name="${base%_trimmed}"
    R2=""
  else
    R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
      continue
    fi
    sample_name="$(basename "$R1" _R1_val_1.fq.gz)"
  fi

  marked_bam="${BAM_DIR}/${sample_name}_marked.bam"
  sorted_bam="${BAM_DIR}/${sample_name}_sorted.bam"
  shifted_bam="${BAM_DIR}/${sample_name}_sorted.shifted.bam"
  bw_file="${TRACK_DIR}/${sample_name}.bw"
  bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

  # 1) Align + dedup -> _marked.bam (keep dup-marked BAM for preseq)
  if [[ -f "$marked_bam" && -f "${marked_bam}.bai" ]]; then
    echo "[STEP2] Found ${sample_name}_marked.bam, skipping alignment+markdup."
  else
    echo "[STEP2] Aligning + markdup: $sample_name"
    align_and_dedup "$R1" "$R2" "$marked_bam" "$bt2_log"
  fi

  # 2) Filter -> _sorted.bam (dups + chrM + low-MAPQ removed)
  if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
    echo "[STEP2] Found ${sample_name}_sorted.bam, skipping filter."
  else
    echo "[STEP2] Filter (MAPQ >= $MAPQ_MIN, exclude $EXCLUDE_CHROMS, drop dups): $sample_name"
    filter_bam "$marked_bam" "$sorted_bam"
  fi

  # 3) BAM QC stats from filtered BAM (matches what downstream uses)
  [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
  [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"
  # Pre-filter idxstats too — useful for chrM% / mapping-rate diagnostics
  [[ ! -f "${BAMQC_DIR}/${sample_name}.marked.idxstats.txt" ]] && samtools idxstats "$marked_bam" > "${BAMQC_DIR}/${sample_name}.marked.idxstats.txt"

  # 4) Tn5 shift
  if [[ "$TN5_SHIFT" -eq 1 ]]; then
    if [[ -f "$shifted_bam" && -f "${shifted_bam}.bai" ]]; then
      echo "[STEP2] Found ${sample_name}_sorted.shifted.bam, skipping Tn5 shift."
    else
      echo "[STEP2] Tn5 shift (alignmentSieve --ATACshift): $sample_name"
      tn5_shift "$sorted_bam" "$shifted_bam"
    fi
    track_bam="$shifted_bam"
  else
    track_bam="$sorted_bam"
  fi

  # 5) bigWig from shifted BAM
  if [[ ! -f "$bw_file" ]]; then
    bam_cov_args=(-b "$track_bam" -o "$bw_file" -p "$BAMCOV_CPU" -bs "$BINSIZE" \
      --effectiveGenomeSize "$GENOMESIZE" --normalizeUsing "$NORMALIZE" \
      --ignoreForNormalization "$IGNORE_CHR")
    [[ -f "${BLACKLIST:-}" ]] && bam_cov_args+=( --blackListFileName "$BLACKLIST" )
    bamCoverage "${bam_cov_args[@]}"
  fi

  # 6) Optional cleanup of marked BAM (preseq won't have run yet on first pass)
  if [[ "$DELETE_MARKED_AFTER_BOWTIE2" -eq 1 ]]; then
    rm -f "$marked_bam" "${marked_bam}.bai"
  fi
done

if [[ "$DELETE_TRIMMED_AFTER_ALIGN" -eq 1 ]]; then
  if [[ "${SE:-0}" -eq 1 ]]; then
    rm -f "${TRIM_DIR}"/*_trimmed.fq.gz || true
  else
    rm -f "${TRIM_DIR}"/*_val_*.fq.gz || true
  fi
fi

echo "=== STEP 2 complete ==="
