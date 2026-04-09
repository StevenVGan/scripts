#!/usr/bin/env bash
# Symlink raw FASTQs into project data/ from a TSV map (Illumina *_R1/_R2_001 or SE SRR.fastq.gz).
# Location: pipeline/tools/prep/link_fastq.sh
set -euo pipefail

# ==== CONFIG: Edit defaults below, or export RAW_DIR / DEST_DIR / MAP_FILE =====
# (project wrappers can export those vars and exec this script unchanged)
RAW_DIR="${RAW_DIR:-$HOME/work/raw_seq/YOUR_RUN_ID}"
DEST_DIR="${DEST_DIR:-$HOME/work/seq/CUTRUN/YOUR_PROJECT/data}"
MAP_FILE="${MAP_FILE:-$HOME/work/seq/CUTRUN/YOUR_PROJECT/link_sample.tsv}"
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

  # pattern: SG388*_R1_001.fastq.gz (supports both SG388_S1_L001_... and SG388_CnR_..._S1_L001_...)
  # SE fallback: ENA/SRA-style ${prefix}.fastq.gz (e.g. SRR123.fastq.gz) when Illumina R1 name not found
  pattern_R1="${RAW_DIR}/${prefix}"*"_R1_001.fastq.gz"
  pattern_R2="${RAW_DIR}/${prefix}"*"_R2_001.fastq.gz"
  se_candidate="${RAW_DIR}/${prefix}.fastq.gz"

  # nullglob: if no Illumina-named files exist, do not treat the glob string as one path
  shopt -s nullglob
  matches_R1=( $pattern_R1 )
  matches_R2=( $pattern_R2 )
  shopt -u nullglob

  if (( ${#matches_R1[@]} == 1 )); then
    src_R1="${matches_R1[0]}"
  elif [[ -f "$se_candidate" ]]; then
    src_R1="$se_candidate"
  else
    echo "WARNING: ${prefix} R1: expected 1 Illumina *_R1_001 match or ${prefix}.fastq.gz, got ${#matches_R1[@]} (pattern: $pattern_R1)" >&2
  fi

  if (( ${#matches_R2[@]} == 1 )); then
    src_R2="${matches_R2[0]}"
  else
    echo "WARNING: ${prefix} R2: expected 1 match, got ${#matches_R2[@]} (pattern: $pattern_R2)" >&2
  fi

  # ---- Determine destination paths ------------------------------------------
  dest_R1="${DEST_DIR}/${newname}_R1.fastq.gz"
  dest_R2="${DEST_DIR}/${newname}_R2.fastq.gz"

  # ---- Create symlinks if sources exist (force replace dest) ---------------
  if [[ -n "$src_R1" ]]; then
    echo "ln -sfn \"$src_R1\" \"$dest_R1\""
    ln -sfn "$src_R1" "$dest_R1"
  fi

  if [[ -n "$src_R2" ]]; then
    echo "ln -sfn \"$src_R2\" \"$dest_R2\""
    ln -sfn "$src_R2" "$dest_R2"
  else
    rm -f "$dest_R2"
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
