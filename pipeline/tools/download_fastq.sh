#!/usr/bin/env bash
# Download CnR/Seq FASTQ data from IGM FTP (e.g. submission 1027).
# Usage: export FTP_PASSWORD='your_password'; ./download_fastq.sh
set -euo pipefail

# ==== CONFIG: Edit these for your project =====================================
# Example: CnR Seq data from submission 1027
# Sequencing data URL (FTP base, no trailing slash)
FTP_BASE="ftp://igm-storage.ucsd.edu/260212_LH00444_0470_A233NMNLT3"
# FTP credentials (password often contains md5sum from sequencing center)
FTP_USER="rosenfeld"
# Set via env for security: export FTP_PASSWORD='2Zlgavv9vj'
# Or set here (not recommended for shared repos):
# FTP_PASSWORD="${FTP_PASSWORD:-}"
FTP_PASSWORD="2Zlgavv9vj"

# Destination: raw_seq for downloaded data, or project-specific
DEST_DIR="$HOME/work/raw_seq/260212_LH00444_0470_A233NMNLT3"
# File prefix filter: sequencing data typically start with SG
FILE_PREFIX="${FILE_PREFIX:-SG}"
# ==============================================================================

if [[ -z "$FTP_PASSWORD" ]]; then
  echo "ERROR: Set FTP_PASSWORD (e.g. export FTP_PASSWORD='your_password')" >&2
  exit 1
fi

mkdir -p "$DEST_DIR"

wget -r -nH -nd -P "$DEST_DIR" --no-passive-ftp \
  -A "${FILE_PREFIX}*,*md5*" -R "index.html*" \
  --user="$FTP_USER" --password="$FTP_PASSWORD" \
  "${FTP_BASE}/"

# Verify integrity if md5 checksum file(s) were downloaded
for md5file in "$DEST_DIR"/*md5*; do
  [[ -f "$md5file" ]] || continue
  echo "Verifying: $(basename "$md5file")" >&2
  (cd "$DEST_DIR" && md5sum -c --ignore-missing "$(basename "$md5file")")
done

echo "Download complete." >&2
