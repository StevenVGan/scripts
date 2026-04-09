#!/usr/bin/env bash
# Merge Illumina per-lane PE FASTQs in place, then delete the lane files.
#
# Required env (or export before running):
#   IGM_DIR     — directory containing lane FASTQs (e.g. raw_seq/RUN_ID)
#   MAP_FILE    — TSV: col1 = basename prefix matching ${prefix}_*_R1_001.fastq.gz, col2 = output sample name
#
# Outputs: ${IGM_DIR}/${sample_name}_R1.fastq.gz and _R2.fastq.gz
#
# Example:
#   IGM_DIR=~/work/raw_seq/260401_LH00444_0498_B235NNMLT4 \
#   MAP_FILE=~/work/seq/CUTRUN/my_project/igm_to_sample.tsv \
#     ./merge_lanes_inplace.sh
set -euo pipefail

: "${IGM_DIR:?Set IGM_DIR to the directory with per-lane FASTQs}"
: "${MAP_FILE:?Set MAP_FILE to the igm_prefix<TAB>sample_name TSV}"

if [[ ! -d "$IGM_DIR" ]]; then
  echo "ERROR: IGM_DIR not a directory: $IGM_DIR" >&2
  exit 1
fi
if [[ ! -f "$MAP_FILE" ]]; then
  echo "ERROR: MAP_FILE not found: $MAP_FILE" >&2
  exit 1
fi

while read -r igm_prefix sample_name _rest; do
  [[ -z "${igm_prefix:-}" ]] && continue
  [[ "$igm_prefix" =~ ^# ]] && continue

  shopt -s nullglob
  r1=( "${IGM_DIR}/${igm_prefix}"_*_R1_001.fastq.gz )
  r2=( "${IGM_DIR}/${igm_prefix}"_*_R2_001.fastq.gz )
  shopt -u nullglob

  if (( ${#r1[@]} == 0 )); then
    echo "ERROR: no R1 files matching ${IGM_DIR}/${igm_prefix}_*_R1_001.fastq.gz" >&2
    exit 1
  fi
  if (( ${#r2[@]} == 0 )); then
    echo "ERROR: no R2 files matching ${IGM_DIR}/${igm_prefix}_*_R2_001.fastq.gz" >&2
    exit 1
  fi

  mapfile -t r1s < <(printf '%s\n' "${r1[@]}" | sort)
  mapfile -t r2s < <(printf '%s\n' "${r2[@]}" | sort)

  out_r1="${IGM_DIR}/${sample_name}_R1.fastq.gz"
  out_r2="${IGM_DIR}/${sample_name}_R2.fastq.gz"

  echo "[merge] ${sample_name}: cat ${#r1s[@]} R1 -> $(basename "$out_r1")" >&2
  cat "${r1s[@]}" > "$out_r1"
  echo "[merge] ${sample_name}: cat ${#r2s[@]} R2 -> $(basename "$out_r2")" >&2
  cat "${r2s[@]}" > "$out_r2"

  gzip -t "$out_r1" "$out_r2"

  echo "[merge] ${sample_name}: removing ${#r1s[@]} R1 + ${#r2s[@]} R2 lane files" >&2
  rm -f -- "${r1s[@]}" "${r2s[@]}"
done < "$MAP_FILE"

shopt -s nullglob
leftover=( "${IGM_DIR}"/*_L0*_R1_001.fastq.gz "${IGM_DIR}"/*_L0*_R2_001.fastq.gz )
shopt -u nullglob
if (( ${#leftover[@]} > 0 )); then
  echo "WARN: lane-style FASTQs still present (not in map or pattern mismatch?):" >&2
  printf '  %s\n' "${leftover[@]}" >&2
fi

echo "Done. Merged FASTQs in: ${IGM_DIR}" >&2
