#!/usr/bin/env bash
# Download FASTQs for SRA run accessions (SRR*) from ENA over HTTPS (curl).
# Typical source: GEO/SRA metadata compiled into a one-SRR-per-line list.
# Use when NCBI prefetch/SRA toolkit fails (e.g. TLS issues on some clusters).
#
# SE-first: ENA fastq_ftp is assumed to be a single file per SRR. Paired-end
# libraries may expose two URLs (comma-separated); this script does not split them yet.
#
# Parallel transfers: DOWNLOAD_JOBS (default 8), e.g. DOWNLOAD_JOBS=12 ./download_geo_fastq_ena.sh
#
# Overrides: DEST_DIR=... SRR_LIST_FILE=... LOG_DIR=... MD5_FILE=... ./download_geo_fastq_ena.sh
# Optional positional: ./download_geo_fastq_ena.sh /path/to/srr_accessions.txt
set -euo pipefail

# ==== CONFIG: edit paths for your project (or set env vars instead) =================
# Where SRR123.fastq.gz files are written
DEST_DIR="${DEST_DIR:-${HOME}/work/raw_seq/MY_GEO_RUN/fastq}"
# One SRR accession per line; # comments and blank lines allowed
SRR_LIST_FILE="${SRR_LIST_FILE:-${HOME}/work/raw_seq/MY_GEO_RUN/srr_accessions.txt}"
# Logs and transient manifest (recreated each run)
LOG_DIR="${LOG_DIR:-$(dirname "$DEST_DIR")/logs}"
# Aggregated md5 lines for md5sum -c (basenames relative to DEST_DIR)
MD5_FILE="${MD5_FILE:-$(dirname "$DEST_DIR")/md5sum_ena.txt}"
# Concurrent curl transfers
DOWNLOAD_JOBS="${DOWNLOAD_JOBS:-8}"
# ==================================================================================

if [[ -n "${1:-}" && -f "$1" ]]; then
  SRR_LIST_FILE="$1"
fi

FQ_DIR="$DEST_DIR"
MD5_PART_DIR="${LOG_DIR}/md5_parts"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
MANIFEST="${LOG_DIR}/ena_manifest.tsv"

if [[ ! -f "$SRR_LIST_FILE" ]]; then
  echo "ERROR: SRR list not found: $SRR_LIST_FILE" >&2
  echo "       Edit CONFIG in this script or set SRR_LIST_FILE=... (or pass path as first argument)." >&2
  exit 1
fi

mkdir -p "$FQ_DIR" "$LOG_DIR" "$MD5_PART_DIR"
rm -f "${MD5_PART_DIR}"/*.txt 2>/dev/null || true
: > "${MD5_FILE}.new"

mapfile -t SRRS < <(grep -v '^#' "$SRR_LIST_FILE" | grep -v '^[[:space:]]*$' || true)
if (( ${#SRRS[@]} == 0 )); then
  echo "ERROR: No SRRs in $SRR_LIST_FILE" >&2
  exit 1
fi

: > "$MANIFEST"

echo "[1/2] Resolving ENA URLs (sequential, small requests)..."
for srr in "${SRRS[@]}"; do
  out="${FQ_DIR}/${srr}.fastq.gz"
  if [[ -f "$out" && -s "$out" ]]; then
    echo "[SKIP resolve] $srr (already have FASTQ)"
    continue
  fi
  line="$(curl -sS "${ENA_API}?accession=${srr}&fields=fastq_ftp,fastq_md5&result=read_run&format=tsv" | tail -1)"
  ftp_path="$(echo "$line" | cut -f2)"
  md5="$(echo "$line" | cut -f3)"
  if [[ -z "$ftp_path" || "$ftp_path" == fastq_ftp ]]; then
    echo "ERROR: No fastq_ftp for $srr (line=$line)" >&2
    exit 1
  fi
  if [[ "$ftp_path" == *","* ]]; then
    echo "ERROR: $srr has multiple fastq_ftp entries (likely PE). This tool is SE-only for now." >&2
    exit 1
  fi
  if [[ "$ftp_path" == ftp.sra.ebi.ac.uk/* ]]; then
    url="https://ftp.sra.ebi.ac.uk/${ftp_path#ftp.sra.ebi.ac.uk/}"
  else
    url="https://ftp.sra.ebi.ac.uk/${ftp_path}"
  fi
  printf '%s\t%s\t%s\n' "$srr" "$url" "$md5" >> "$MANIFEST"
done

if [[ ! -s "$MANIFEST" ]]; then
  echo "Nothing to download (all FASTQs present or empty manifest)."
else
  echo "[2/2] Downloading with up to ${DOWNLOAD_JOBS} parallel curl jobs..."
  set -m
  pids=()
  while IFS=$'\t' read -r srr url md5; do
    [[ -z "${srr:-}" ]] && continue
    while (( $(jobs -r | wc -l) >= DOWNLOAD_JOBS )); do
      sleep 0.3
    done
    out="${FQ_DIR}/${srr}.fastq.gz"
    (
      echo "[curl] $srr"
      tmp="${out}.part.${BASHPID}.${RANDOM}"
      rm -f "$tmp"
      if curl -fsSL --connect-timeout 30 --retry 5 --retry-delay 10 \
          -o "$tmp" "$url"; then
        mv -f "$tmp" "$out"
        if [[ -n "$md5" && "$md5" != fastq_md5 ]]; then
          echo "$md5  $(basename "$out")" > "${MD5_PART_DIR}/${srr}.txt"
        fi
      else
        rm -f "$tmp" "$out"
        exit 1
      fi
    ) &
    pids+=( "$!" )
  done < "$MANIFEST"
  dl_fail=0
  for pid in "${pids[@]}"; do
    wait "$pid" || dl_fail=1
  done
  set +m
  if [[ "$dl_fail" -ne 0 ]]; then
    echo "ERROR: One or more downloads failed." >&2
    exit 1
  fi
fi

shopt -s nullglob
for f in "${MD5_PART_DIR}"/*.txt; do
  cat "$f" >> "${MD5_FILE}.new"
done

if [[ -s "${MD5_FILE}.new" ]]; then
  if (cd "$FQ_DIR" && md5sum -c --ignore-missing "${MD5_FILE}.new"); then
    mv "${MD5_FILE}.new" "$MD5_FILE"
  else
    echo "[WARN] md5sum check had issues; see ${MD5_FILE}.new" >&2
  fi
fi

echo "Done. FASTQs in ${FQ_DIR}/ ($(ls -1 "$FQ_DIR"/*.fastq.gz 2>/dev/null | wc -l) files)"
