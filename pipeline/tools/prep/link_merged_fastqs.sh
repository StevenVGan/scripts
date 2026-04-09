#!/usr/bin/env bash
# Symlink merged PE FASTQs (${sample}_R1.fastq.gz / _R2.fastq.gz) into a project data directory.
# Use after merge_lanes_inplace.sh (or any workflow that produces that naming).
#
# Required env:
#   MERGED_FASTQ_DIR  — directory containing ${sample_name}_R{1,2}.fastq.gz
#                       (alias: IGM_FASTQ_DIR if MERGED_FASTQ_DIR unset)
#   DEST_DIR          — e.g. project data/ for 1_trim_qc input
#   MAP_FILE          — same TSV as merge (col1 ignored except comments; col2 = sample_name)
#
# Log: ${MAP_FILE%.tsv}.link.log.tsv
set -euo pipefail

MERGED_FASTQ_DIR="${MERGED_FASTQ_DIR:-${IGM_FASTQ_DIR:-}}"
: "${MERGED_FASTQ_DIR:?Set MERGED_FASTQ_DIR or IGM_FASTQ_DIR to the folder with merged *_R1.fastq.gz}"
: "${DEST_DIR:?Set DEST_DIR to the project raw FASTQ dir (e.g. BASE/data)}"
: "${MAP_FILE:?Set MAP_FILE to the igm_prefix<TAB>sample_name TSV}"

if [[ ! -f "$MAP_FILE" ]]; then
  echo "ERROR: MAP_FILE not found: $MAP_FILE" >&2
  exit 1
fi

mkdir -p "$DEST_DIR"

LOG_FILE="${MAP_FILE%.tsv}.link.log.tsv"
TMP_LOG="${LOG_FILE}.tmp"
RUN_DATE="$(date -Iseconds)"
printf "sample_name\tsrc_R1\tsrc_R2\tdest_R1\tdest_R2\tlinked_at\n" > "$TMP_LOG"

while read -r igm_prefix sample_name _rest; do
  [[ -z "${igm_prefix:-}" ]] && continue
  [[ "$igm_prefix" =~ ^# ]] && continue

  src_R1="${MERGED_FASTQ_DIR}/${sample_name}_R1.fastq.gz"
  src_R2="${MERGED_FASTQ_DIR}/${sample_name}_R2.fastq.gz"
  dest_R1="${DEST_DIR}/${sample_name}_R1.fastq.gz"
  dest_R2="${DEST_DIR}/${sample_name}_R2.fastq.gz"

  if [[ ! -f "$src_R1" ]]; then
    echo "ERROR: missing merged R1: $src_R1" >&2
    exit 1
  fi
  if [[ ! -f "$src_R2" ]]; then
    echo "ERROR: missing merged R2: $src_R2" >&2
    exit 1
  fi

  echo "ln -sfn \"$src_R1\" \"$dest_R1\"" >&2
  ln -sfn "$src_R1" "$dest_R1"
  echo "ln -sfn \"$src_R2\" \"$dest_R2\"" >&2
  ln -sfn "$src_R2" "$dest_R2"

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample_name" "$src_R1" "$src_R2" "$dest_R1" "$dest_R2" "$RUN_DATE" \
    >> "$TMP_LOG"
done < "$MAP_FILE"

mv "$TMP_LOG" "$LOG_FILE"
echo "Done. Symlinks in: ${DEST_DIR}  (log: ${LOG_FILE})" >&2
