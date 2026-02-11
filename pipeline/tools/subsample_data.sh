#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Subsample oridata/*.fastq.gz -> data/: 1M reads per file, pairing preserved.
# Uses seqtk sample (same seed for R1 and R2). Requires: seqtk (e.g. conda install -c bioconda seqtk).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Config: CONFIG_FILE env, or ../cutrun/0_config.sh (when run from pipeline/tools), or same dir
if [[ -n "${CONFIG_FILE:-}" && -f "${CONFIG_FILE}" ]]; then
  source "${CONFIG_FILE}"
elif [[ -f "${SCRIPT_DIR}/../cutrun/0_config.sh" ]]; then
  source "${SCRIPT_DIR}/../cutrun/0_config.sh"
elif [[ -f "${SCRIPT_DIR}/0_config.sh" ]]; then
  source "${SCRIPT_DIR}/0_config.sh"
fi
[[ -n "${BASE:-}" ]] || { echo "ERROR: BASE not set. Set CONFIG_FILE=/path/to/0_config.sh or run from project with config." >&2; exit 1; }

if ! command -v seqtk &>/dev/null; then
  echo "ERROR: seqtk not found. Install with: conda install -c bioconda seqtk" >&2
  exit 1
fi

N=1000000
SEED=42
ORIDATA="${BASE}/oridata"
DATA="${BASE}/data"

for R1 in "${ORIDATA}"/*_R1.fastq.gz; do
  R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
  [[ -f "$R2" ]] || { echo "[WARN] No pair: $R2"; continue; }
  name="$(basename "$R1" _R1.fastq.gz)"
  echo "$name: sampling $N reads (seed=$SEED) -> data/"
  seqtk sample -s "$SEED" "$R1" "$N" | gzip > "${DATA}/${name}_R1.fastq.gz"
  seqtk sample -s "$SEED" "$R2" "$N" | gzip > "${DATA}/${name}_R2.fastq.gz"
done

echo "Done."
