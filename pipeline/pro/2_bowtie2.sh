#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 2_bowtie2.sh — PRO-seq alignment + strand-aware bigWig generation
#
# Forked from csRNA/2_bowtie2.sh. The only change is in run_bamcoverage_tracks():
# when PROSEQ_FLIP_STRAND=1, the --filterRNAstrand labels are SWAPPED so that
# ${sample}_fwd.bw represents the nascent RNA on the + strand (PRO-seq R1 is
# the reverse complement of the nascent RNA; see 0_config.sh notes).
#
# Inputs
#   - Trimmed FASTQs in ${TRIM_DIR} (PE: *_R1_val_1.fq.gz + *_R2_val_2.fq.gz; SE: *_trimmed.fq.gz)
#   - Bowtie2 genome index: ${GENOME_INDEX}
#
# Outputs
#   - ${BAM_DIR}/*_sorted.bam + .bai, *_bowtie2.log
#   - ${BAMQC_DIR}/*.idxstats.txt, *.stats
#   - ${TRACK_DIR}/${sample}.bw (optional unstranded, when COMBINED_BIGWIG=1)
#   - ${TRACK_DIR}/${sample}_fwd.bw and _rev.bw (when STRAND_BIGWIG=1; RNA-strand when PROSEQ_FLIP_STRAND=1)
###############################################################################


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "2_bowtie2"

check_cmd bowtie2
check_cmd samtools
check_cmd bamCoverage

echo "=== STEP 2: Bowtie2 -> sorted BAM + bigWig + basic BAM QC ==="
echo "[INFO] GENOME_INDEX:       $GENOME_INDEX"
echo "[INFO] TRIM_DIR:           $TRIM_DIR"
echo "[INFO] BAM_DIR:            $BAM_DIR"
echo "[INFO] TRACK_DIR:          $TRACK_DIR"
echo "[INFO] Mode:               $([ "${SE:-0}" -eq 1 ] && echo "single-end" || echo "paired-end")"
echo "[INFO] BT2_EXTRA:          ${BT2_EXTRA:-}"
echo "[INFO] STRAND_BIGWIG:      ${STRAND_BIGWIG:-0}  COMBINED_BIGWIG: ${COMBINED_BIGWIG:-1}"
echo "[INFO] PROSEQ_FLIP_STRAND: ${PROSEQ_FLIP_STRAND:-0} (swaps fwd/rev bigWig labels so they reflect RNA strand)"

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

# bigWig(s): optional unstranded + strand-separated.
# When PROSEQ_FLIP_STRAND=1, we swap --filterRNAstrand labels so the output
# file named "_fwd.bw" contains reads from the RNA's + strand (which, for
# PRO-seq, are the reads that aligned to the - strand).
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
  if [[ ( "${COMBINED_BIGWIG:-1}" -ne 1 || -f "$bw_file" ) && -f "$bw_fwd" && -f "$bw_rev" ]]; then
    return                                   # every requested track exists
  fi
  local sf
  sf=$(unstranded_scale_factor "$sample_name" "$sorted_bam") || exit 1
  # strand tracks: same filters/blacklist as `common`, but ONE shared (all-reads) scale factor
  local stranded=( "${common[@]}" )
  stranded[${#stranded[@]}]=--scaleFactor; stranded[${#stranded[@]}]="$sf"
  local i; for i in "${!stranded[@]}"; do
    [[ "${stranded[$i]}" == "--normalizeUsing" ]] && stranded[$((i + 1))]=None
  done

  local strand_for_fwd strand_for_rev
  if [[ "${PROSEQ_FLIP_STRAND:-0}" -eq 1 ]]; then
    # PRO-seq: R1 is antisense of nascent RNA → swap filterRNAstrand labels
    strand_for_fwd="reverse"
    strand_for_rev="forward"
  else
    strand_for_fwd="forward"
    strand_for_rev="reverse"
  fi

  if [[ ! -f "$bw_fwd" ]]; then
    bamCoverage "${stranded[@]}" --filterRNAstrand "$strand_for_fwd" -o "$bw_fwd"
  fi
  if [[ ! -f "$bw_rev" ]]; then
    bamCoverage "${stranded[@]}" --filterRNAstrand "$strand_for_rev" -o "$bw_rev"
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
      bt2_cmd=( bowtie2 -x "$GENOME_INDEX" -U "$R1" --very-sensitive-local -p "$BT2_CPU" )
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
      bt2_cmd=( bowtie2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" --very-sensitive-local -p "$BT2_CPU" )
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
