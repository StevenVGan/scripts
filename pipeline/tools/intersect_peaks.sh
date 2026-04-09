#!/usr/bin/env bash
# Thin wrapper: same as peak_ops.sh --mode intersect (BED + UpSet/Venn via --viz).
# Keeps the short "intersect_peaks.sh OUT.bed A B ..." CLI without typing --mode.
#
# Usage:
#   ./intersect_peaks.sh [OPTIONS] OUTPUT.bed INPUT1 INPUT2 [INPUT3 ...]
#
# Options (passed through to peak_ops.sh):
#   --viz none|upset|venn|both   Default upset (see peak_ops.sh).
#   --slop BP, --genome-sizes F, --merge, --names "N1,N2,…"
#
# For union / distinct / full control, call peak_ops.sh directly.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/peak_ops.sh" --mode intersect "$@"
