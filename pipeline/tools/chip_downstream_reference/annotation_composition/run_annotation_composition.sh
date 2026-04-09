#!/usr/bin/env bash
# Run method2_annotation_composition.py from this directory (conda bio for matplotlib if needed).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 not found" >&2
  exit 1
fi
if ! python3 -c "import matplotlib" 2>/dev/null; then
  # shellcheck source=/dev/null
  source "${HOME}/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate bio 2>/dev/null || true
fi
exec python3 "${SCRIPT_DIR}/method2_annotation_composition.py" \
  --cobinding-root "${SCRIPT_DIR}/.." \
  --chip-root "$(cd "${SCRIPT_DIR}/../../.." && pwd)" \
  --out-dir "${SCRIPT_DIR}/output" \
  "$@"
