#!/usr/bin/env bash
# LEGACY — IGM moved deliveries to Globus during 2026 (last FTP pull here 2026-08-11; port 21 was closed by
# 2026-09-23). Use download_igm_fastq_globus.sh for new runs. Kept for reference.
# Download FASTQ data from IGM FTP (wget + md5sum -c).
# Location: pipeline/tools/prep/download_fastq.sh
# Usage: export FTP_USER='lab_account' FTP_PASSWORD='your_password' FTP_BASE='ftp://igm-storage.ucsd.edu/RUN_ID' \
#        DEST_DIR=$HOME/work/raw_seq/RUN_ID FILE_PREFIX=SG; ./download_fastq.sh
set -euo pipefail

# ==== CONFIG: set via env vars (nothing site-specific is hardcoded in this PUBLIC repo) =====
# Sequencing data URL (FTP base, no trailing slash)
FTP_BASE="${FTP_BASE:-ftp://igm-storage.ucsd.edu/YOUR_RUN_ID}"
# FTP credentials — both env-supplied (the account name and the password are the lab's)
FTP_USER="${FTP_USER:-}"
# Set via env for security (never hardcode in this shared/public repo):
#   export FTP_PASSWORD='your_password'
FTP_PASSWORD="${FTP_PASSWORD:-}"

# Destination: raw_seq for downloaded data, or project-specific
DEST_DIR="${DEST_DIR:-$HOME/work/raw_seq/YOUR_RUN_ID}"
# Some IGM CnR runs name FASTQs PS###_... ; other runs often use SG###
FILE_PREFIX="${FILE_PREFIX:-PS}"
# ==============================================================================

if [[ -z "$FTP_USER" || -z "$FTP_PASSWORD" ]]; then
  echo "ERROR: Set FTP_USER and FTP_PASSWORD (e.g. export FTP_USER='lab_account' FTP_PASSWORD='your_password')" >&2
  exit 1
fi

mkdir -p "$DEST_DIR"

wget -r -nH -nd -P "$DEST_DIR" --no-passive-ftp \
  -A "${FILE_PREFIX}*,*md5*" -R "index.html*" \
  --user="$FTP_USER" --password="$FTP_PASSWORD" \
  "${FTP_BASE}/"

# Verify integrity if md5 checksum file(s) were downloaded. Every listed file must be present (checked explicitly —
# coreutils 8.25's `md5sum --ignore-missing` can silently drop an existing file), then plain md5sum -c.
for md5file in "$DEST_DIR"/*md5*; do
  [[ -f "$md5file" ]] || continue
  echo "Verifying: $(basename "$md5file")" >&2
  missing=$(grep -vE '^\s*(#|$)' "$md5file" | sed -E 's/^[0-9A-Fa-f]{32}[ *]+//' | while IFS= read -r f; do [[ -e "$DEST_DIR/$f" ]] || printf '%s\n' "$f"; done)
  if [[ -n "$missing" ]]; then echo "ERROR: listed in $(basename "$md5file") but not downloaded:" >&2; printf '  %s\n' $missing >&2; exit 1; fi
  (cd "$DEST_DIR" && md5sum -c "$(basename "$md5file")")
done

echo "Download complete." >&2
