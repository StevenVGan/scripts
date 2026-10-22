#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 2_bowtie2.sh
#
# Purpose
#   Step 2 of the pipeline: alignment and BAM/track generation.
#   Performs Bowtie2 alignment, produces sorted+indexed BAMs, generates bigWig
#   tracks, and writes basic BAM QC stats (idxstats, samtools stats).
#
# Inputs
#   - Trimmed paired FASTQs in ${TRIM_DIR}:
#       *_R1_val_1.fq.gz and *_R2_val_2.fq.gz
#   - Bowtie2 genome index: ${GENOME_INDEX} (configured in 0_config.sh)
#
# Outputs
#   - Sorted BAM + index:
#       ${BAM_DIR}/*_sorted.bam and ${BAM_DIR}/*_sorted.bam.bai
#   - Bowtie2 logs:
#       ${BAM_DIR}/*_bowtie2.log
#   - Per-contig mapping stats:
#       ${BAMQC_DIR}/*.idxstats.txt
#   - samtools stats (for MultiQC):
#       ${BAMQC_DIR}/*.stats
#   - bigWig tracks (see STRAND_BIGWIG / COMBINED_BIGWIG in 0_config.sh):
#       ${TRACK_DIR}/${sample}.bw (optional unstranded)
#       ${TRACK_DIR}/${sample}_fwd.bw and _rev.bw (strand-separated when STRAND_BIGWIG=1)
#
# Behavior
#   - Streams Bowtie2 output directly into samtools sort (does not write SAM).
#   - Skips alignment if *_sorted.bam + .bai already exist.
#   - Skips bamCoverage if the .bw already exists.
#   - Optional: deletes trimmed FASTQs after alignment (space-saving), controlled by
#       DELETE_TRIMMED_AFTER_ALIGN in 0_config.sh.
#
# Logging
#   - Writes a timestamped log to ${LOG_DIR}/2_bowtie2_*.log
#
# Requirements
#   - bowtie2, samtools, bamCoverage on PATH.
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "2_bowtie2"

check_cmd bowtie2
check_cmd samtools
check_cmd bamCoverage

echo "=== STEP 2: Bowtie2 -> sorted BAM + bigWig + basic BAM QC ==="
echo "[INFO] GENOME_INDEX: $GENOME_INDEX"
echo "[INFO] TRIM_DIR:     $TRIM_DIR"
echo "[INFO] BAM_DIR:      $BAM_DIR"
echo "[INFO] TRACK_DIR:   $TRACK_DIR"
echo "[INFO] Mode:        $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"
echo "[INFO] BT2_PRESET:  ${BT2_PRESET:---very-sensitive-local}"
echo "[INFO] BT2_EXTRA:   ${BT2_EXTRA:-}"
echo "[INFO] STRAND_BIGWIG=${STRAND_BIGWIG:-0} COMBINED_BIGWIG=${COMBINED_BIGWIG:-1}"

# Scale factor of the UNSTRANDED (all-reads) track, from deepTools' own calculation (bamCoverage --verbose
# prints "Final scaling factor"). Both strand files reuse it: with --normalizeUsing, deepTools computes the
# factor on each strand's OWN filtered reads, which puts + and − on different scales (csRNA 260915: − file ~2×
# the + file, because structural RNA puts ~2/3 of reads on one strand). Cached in <sample>.scale_factor, reused
# only while newer than the BAM and the combined track is not being rebuilt. Reads `common` from the caller (bash dynamic scope). Prints the factor on stdout.
unstranded_scale_factor() {
  local sample_name="$1"
  local sorted_bam="$2"
  local sf_file="${TRACK_DIR}/${sample_name}.scale_factor"
  local bw_file="${TRACK_DIR}/${sample_name}.bw"
  local bw_fwd="${TRACK_DIR}/${sample_name}_fwd.bw" bw_rev="${TRACK_DIR}/${sample_name}_rev.bw"
  local want_combined=0
  [[ "${COMBINED_BIGWIG:-1}" -eq 1 && ! -f "$bw_file" ]] && want_combined=1
  if [[ -z "${NORMALIZE:-}" || "$NORMALIZE" == "None" ]]; then
    # no normalisation: the factor is 1, but the combined track is still built when it is wanted
    if (( want_combined )); then
      bamCoverage "${common[@]}" -o "$bw_file" >&2 || { echo "[ERROR] bamCoverage failed for ${sample_name}" >&2; rm -f "$bw_file"; return 1; }
    fi
    echo 1; return
  fi
  if [[ -s "$sf_file" && "$sf_file" -nt "$sorted_bam" ]] && (( ! want_combined )); then
    cat "$sf_file"; return
  fi
  # Cache miss. Build the combined track if it is wanted; otherwise a throwaway probe restricted to the shortest
  # contig >= 10 Mb (deepTools derives the factor from whole-BAM idxstats before any region work, so --region only
  # shortens the coverage pass; the probe lands in TMPDIR, never beside the real tracks). Strand tracks built under
  # an older factor are removed first so a resume can never leave the pair on two scales.
  rm -f "$bw_fwd" "$bw_rev"
  local out="$bw_file" probe=0 region=() contig=""
  if (( ! want_combined )); then
    out="${TMPDIR:-/tmp}/${sample_name}.sf_probe.$$.bw"; probe=1
    contig="$(samtools view -H "$sorted_bam" | awk -F'\t' '$1=="@SQ"{n="";l=0; for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4); if($i~/^LN:/)l=substr($i,4)+0}; if(l>=10000000 && (best==""||l<bl)){best=n;bl=l}} END{print best}')"
    if [[ -n "$contig" ]]; then region=( --region "$contig" ); fi
  fi
  local vlog="${TRACK_DIR}/${sample_name}.bamCoverage_verbose.log"
  if ! bamCoverage "${common[@]}" ${region[@]+"${region[@]}"} --verbose -o "$out" > "$vlog" 2>&1; then
    echo "[ERROR] bamCoverage failed for ${sample_name} — see $vlog" >&2
    rm -f "$out"
    return 1
  fi
  local sf
  sf=$(awk -F': ' '/^Final scaling factor:/{v=$2} END{print v}' "$vlog")
  if [[ -z "$sf" ]]; then
    echo "[ERROR] no 'Final scaling factor' in $vlog" >&2
    rm -f "$out"
    return 1
  fi
  cat "$vlog" >&2        # keep deepTools' own messages (blacklist counts, warnings) in the step log
  rm -f "$vlog"
  (( probe )) && rm -f "$out"
  echo "$sf" > "$sf_file"
  echo "[INFO] ${sample_name}: unstranded ${NORMALIZE} scale factor ${sf} (applied to both strand tracks)" >&2
  echo "$sf"
}

# bigWig(s): optional unstranded + strand-separated (deeptools --filterRNAstrand)
run_bamcoverage_tracks() {
  local sample_name="$1"
  local sorted_bam="$2"
  local common=( -b "$sorted_bam" -p "$BAMCOV_CPU" -bs "$BINSIZE" --effectiveGenomeSize "$GENOMESIZE" --normalizeUsing "$NORMALIZE" --ignoreForNormalization "$IGNORE_CHR" )
  [[ "${BAMCOV_IGNORE_DUP:-0}" -eq 1 ]] && common+=( --ignoreDuplicates )
  [[ -f "${BLACKLIST:-}" ]] && common+=( --blackListFileName "$BLACKLIST" )

  if [[ "${STRAND_BIGWIG:-0}" -ne 1 ]]; then
    local bw_file="${TRACK_DIR}/${sample_name}.bw"
    if [[ ! -f "$bw_file" ]]; then
      bamCoverage "${common[@]}" -o "$bw_file"
    fi
    return
  fi

  local bw_file="${TRACK_DIR}/${sample_name}.bw"
  local bw_fwd="${TRACK_DIR}/${sample_name}_fwd.bw"
  local bw_rev="${TRACK_DIR}/${sample_name}_rev.bw"
  if [[ "${COMBINED_BIGWIG:-1}" -eq 1 && ! -f "$bw_file" ]] || [[ ! -f "$bw_fwd" || ! -f "$bw_rev" ]]; then
    local sf
    sf=$(unstranded_scale_factor "$sample_name" "$sorted_bam") || exit 1
    # strand tracks: same filters/blacklist as `common`, but ONE shared (all-reads) scale factor
    local stranded=( "${common[@]}" )
    stranded[${#stranded[@]}]=--scaleFactor; stranded[${#stranded[@]}]="$sf"
    local i; for i in "${!stranded[@]}"; do
      [[ "${stranded[$i]}" == "--normalizeUsing" ]] && stranded[$((i + 1))]=None
    done
    if [[ ! -f "$bw_fwd" ]]; then
      bamCoverage "${stranded[@]}" --filterRNAstrand forward -o "$bw_fwd"
    fi
    if [[ ! -f "$bw_rev" ]]; then
      bamCoverage "${stranded[@]}" --filterRNAstrand reverse -o "$bw_rev"
    fi
  fi
}

if ! compgen -G "${GENOME_INDEX}.1.bt2*" > /dev/null; then
  echo "[ERROR] bowtie2 index not found at ${GENOME_INDEX}(.1.bt2) — is the NAS mounted / REF_ROOT set for a no-NAS node? (REF_ROOT=${REF_ROOT})" >&2
  exit 1
fi

if [[ "${SE:-0}" -eq 1 ]]; then
  # Single-end: *_trimmed.fq.gz, bowtie2 -U
  trimmed_files=( "${TRIM_DIR}"/*_trimmed.fq.gz )
  if (( ${#trimmed_files[@]} == 0 )); then
    echo "[ERROR] No *_trimmed.fq.gz files found in $TRIM_DIR"
    exit 1
  fi
  for R1 in "${trimmed_files[@]}"; do
    base="$(basename "$R1" .fq.gz)"
    if [[ "$base" == *_R1_trimmed ]]; then
      sample_name="${base%_R1_trimmed}"
    else
      sample_name="${base%_trimmed}"
    fi
    sorted_bam="${BAM_DIR}/${sample_name}_sorted.bam"
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
      echo "[STEP2] Found existing BAM+BAI for $sample_name, skipping alignment."
    else
      echo "[STEP2] Aligning (single-end): $sample_name"
      bt2_cmd=( bowtie2 -x "$GENOME_INDEX" -U "$R1" "${BT2_PRESET:---very-sensitive-local}" -p "$BT2_CPU" )
      if [[ -n "${BT2_EXTRA:-}" ]]; then
        read -r -a bt2_extra_arr <<< "$BT2_EXTRA"
        bt2_cmd+=( "${bt2_extra_arr[@]}" )
      fi
      "${bt2_cmd[@]}" 2> "$bt2_log" \
        | samtools sort -@ "$BT2_CPU" -o "$sorted_bam" -
      samtools index "$sorted_bam"
    fi

    [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
    [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"

    run_bamcoverage_tracks "$sample_name" "$sorted_bam"
  done
  [[ "$DELETE_TRIMMED_AFTER_ALIGN" -eq 1 ]] && rm -f "${TRIM_DIR}"/*_trimmed.fq.gz || true
else
  # Paired-end
  r1_files=( "${TRIM_DIR}"/*_R1_val_1.fq.gz )
  if (( ${#r1_files[@]} == 0 )); then
    echo "[ERROR] No *_R1_val_1.fq.gz files found in $TRIM_DIR"
    exit 1
  fi

  for R1 in "${r1_files[@]}"; do
    R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "[WARN] Missing pair for $R1 (expected $R2). Skipping."
      continue
    fi

    sample_name="$(basename "$R1" _R1_val_1.fq.gz)"
    sorted_bam="${BAM_DIR}/${sample_name}_sorted.bam"
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    if [[ -f "$sorted_bam" && -f "${sorted_bam}.bai" ]]; then
      echo "[STEP2] Found existing BAM+BAI for $sample_name, skipping alignment."
    else
      echo "[STEP2] Aligning: $sample_name"
      bt2_cmd=( bowtie2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" "${BT2_PRESET:---very-sensitive-local}" -p "$BT2_CPU" )
      if [[ -n "${BT2_EXTRA:-}" ]]; then
        read -r -a bt2_extra_arr <<< "$BT2_EXTRA"
        bt2_cmd+=( "${bt2_extra_arr[@]}" )
      fi
      "${bt2_cmd[@]}" 2> "$bt2_log" \
        | samtools sort -@ "$BT2_CPU" -o "$sorted_bam" -
      samtools index "$sorted_bam"
    fi

    [[ ! -f "${BAMQC_DIR}/${sample_name}.idxstats.txt" ]] && samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"
    [[ ! -f "${BAMQC_DIR}/${sample_name}.stats" ]] && samtools stats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.stats"

    run_bamcoverage_tracks "$sample_name" "$sorted_bam"
  done

  [[ "$DELETE_TRIMMED_AFTER_ALIGN" -eq 1 ]] && rm -f "${TRIM_DIR}"/*_val_*.fq.gz || true
fi

echo "=== STEP 2 complete ==="
