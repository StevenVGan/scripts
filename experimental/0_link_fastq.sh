#!/usr/bin/env bash
set -euo pipefail

# ==== CONFIG (adjust for each run) ============================================
RAW_DIR="$HOME/work/raw_seq/260115_LH00444_0462_A233NNNLT3"
DEST_DIR="$HOME/work/seq/CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka/data"
MAP_FILE="$HOME/work/seq/CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka/link_sample.tsv"
# ==============================================================================

mkdir -p "$DEST_DIR"

# New: separate log file instead of overwriting MAP_FILE
LOG_FILE="${MAP_FILE%.tsv}.log.tsv"
TMP_LOG="${LOG_FILE}.tmp"
RUN_DATE="$(date -Iseconds)"

# Write header of the enriched log
printf "prefix\tnewname\tsrc_R1\tsrc_R2\tdest_R1\tdest_R2\tlinked_at\n" > "$TMP_LOG"

while read -r prefix newname _rest; do
  # skip empty lines or comment lines
  [[ -z "${prefix:-}" ]] && continue
  [[ "$prefix" =~ ^# ]] && continue

  # ---- Resolve R1/R2 source paths -------------------------------------------
  src_R1=""
  src_R2=""

  # pattern: SG387_CnR_IgG_MCF7_DMSO_E2_S15_L001_R2_001.fastq.gz
  pattern_R1="${RAW_DIR}/${prefix}"*"R1_001.fastq.gz"
  pattern_R2="${RAW_DIR}/${prefix}"*"R2_001.fastq.gz"

  matches_R1=( $pattern_R1 )
  matches_R2=( $pattern_R2 )

  if (( ${#matches_R1[@]} == 1 )); then
    src_R1="${matches_R1[0]}"
  else
    echo "WARNING: ${prefix} R1: expected 1 match, got ${#matches_R1[@]} (pattern: $pattern_R1)" >&2
  fi

  if (( ${#matches_R2[@]} == 1 )); then
    src_R2="${matches_R2[0]}"
  else
    echo "WARNING: ${prefix} R2: expected 1 match, got ${#matches_R2[@]} (pattern: $pattern_R2)" >&2
  fi

  # ---- Determine destination paths ------------------------------------------
  dest_R1="${DEST_DIR}/${newname}_R1.fastq.gz"
  dest_R2="${DEST_DIR}/${newname}_R2.fastq.gz"

  # ---- Create symlinks if sources exist -------------------------------------
  if [[ -n "$src_R1" ]]; then
    if [[ -e "$dest_R1" ]]; then
      echo "SKIP: $dest_R1 already exists" >&2
    else
      echo "ln -s \"$src_R1\" \"$dest_R1\""
      ln -s "$src_R1" "$dest_R1"
    fi
  fi

  if [[ -n "$src_R2" ]]; then
    if [[ -e "$dest_R2" ]]; then
      echo "SKIP: $dest_R2 already exists" >&2
    else
      echo "ln -s \"$src_R2\" \"$dest_R2\""
      ln -s "$src_R2" "$dest_R2"
    fi
  fi

  # ---- Append enriched line to TMP_LOG --------------------------------------
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$prefix" \
    "$newname" \
    "${src_R1:-"-"}" \
    "${src_R2:-"-"}" \
    "$dest_R1" \
    "$dest_R2" \
    "$RUN_DATE" \
    >> "$TMP_LOG"

done < "$MAP_FILE"

# Atomically finalize the log file (MAP_FILE is untouched)
mv "$TMP_LOG" "$LOG_FILE"
