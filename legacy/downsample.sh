#!/usr/bin/env bash
set -euo pipefail

# Downsample paired-end FASTQ.gz to exactly N read PAIRS per sample, randomly (pair-preserving).
# Input:  ./oridata/*_R1.fastq.gz and matching *_R2.fastq.gz
# Output: ./data/<sample>_R1.fastq.gz and ./data/<sample>_R2.fastq.gz
#
# Requirements:
#   - seqtk (conda/mamba: mamba install -c conda-forge -c bioconda seqtk)
#   - gzip, awk
#
# Notes:
#   - Random sampling is done from R1, then the same read IDs are used to extract from both R1 and R2.
#   - This preserves pairing.
#   - If a sample has fewer than N pairs, the script will error out (change behavior if you prefer).

# TODO: this is not run yet.

IN_DIR="./oridata"
OUT_DIR="./data"
N_READ_PAIRS=5000000          # target read pairs per sample
SEED=100                      # reproducible randomness

mkdir -p "${OUT_DIR}"

shopt -s nullglob

# Loop over all R1 files and infer sample prefix
for r1 in "${IN_DIR}"/*_R1.fastq.gz; do
  base="$(basename "${r1}")"
  sample="${base%_R1.fastq.gz}"
  r2="${IN_DIR}/${sample}_R2.fastq.gz"

  if [[ ! -f "${r2}" ]]; then
    echo "[WARN] Missing R2 for sample '${sample}': ${r2} (skipping)" >&2
    continue
  fi

  out_r1="${OUT_DIR}/${sample}_R1.fastq.gz"
  out_r2="${OUT_DIR}/${sample}_R2.fastq.gz"

  echo "[INFO] Processing sample: ${sample}"

  # Optional sanity check: ensure enough reads exist in R1
  # (FASTQ has 4 lines per read)
  total_pairs=$(( $(zcat "${r1}" | wc -l) / 4 ))
  if (( total_pairs < N_READ_PAIRS )); then
    echo "[ERROR] Sample '${sample}' has only ${total_pairs} read pairs in R1 (< ${N_READ_PAIRS})." >&2
    exit 1
  fi

  # Make a temporary file for IDs to keep (auto-cleaned)
  keep_ids="$(mktemp)"
  trap 'rm -f "${keep_ids}"' EXIT

  # 1) Randomly sample EXACTLY N reads from R1 using seqtk.
  # 2) Extract read IDs from the sampled FASTQ (header lines only: every 4th line starting at 1).
  # 3) Normalize ID:
  #    - Remove leading '@'
  #    - Remove trailing /1 or /2 if present
  #    - Keep only the first whitespace-delimited token (seqtk subseq expects that)
  seqtk sample -s"${SEED}" "${r1}" "${N_READ_PAIRS}" \
    | awk 'NR%4==1{
        id=$1; sub(/^@/,"",id);
        sub(/\/[12]$/,"",id);
        print id
      }' > "${keep_ids}"

  # 2) Subselect BOTH mates by the same ID list (pair-preserving)
  seqtk subseq "${r1}" "${keep_ids}" | gzip > "${out_r1}"
  seqtk subseq "${r2}" "${keep_ids}" | gzip > "${out_r2}"

  # Clean tmp ids for this sample (trap will also clean on exit)
  rm -f "${keep_ids}"
  trap - EXIT

  # Post-check: confirm output read counts match and equal N
  out_n1=$(( $(zcat "${out_r1}" | wc -l) / 4 ))
  out_n2=$(( $(zcat "${out_r2}" | wc -l) / 4 ))
  if (( out_n1 != N_READ_PAIRS || out_n2 != N_READ_PAIRS )); then
    echo "[ERROR] Output counts mismatch for '${sample}': R1=${out_n1}, R2=${out_n2}, expected=${N_READ_PAIRS}" >&2
    exit 1
  fi

  echo "[OK] Wrote: ${out_r1} and ${out_r2} (${N_READ_PAIRS} read pairs each)"
done

echo "[DONE] All samples processed. Output in: ${OUT_DIR}"
