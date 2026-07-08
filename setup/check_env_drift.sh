#!/usr/bin/env bash
# Advisory: warn when a BUILT conda env has changed since its committed
# scripts/env/<env>.yml was last re-exported (the PROVISION §5.2 ritual). Nothing
# else detects this, and a stale yml means the other node rebuilds an outdated
# env -> different tool versions -> divergent results, silently.
#
# Same-node heuristic (deliberately NOT a cross-node `conda env export` diff,
# which is noisy: exports embed build strings + machine-specific prefix paths):
# compare the newest conda-meta/*.json mtime (conda writes one per installed
# package) against the committed yml's mtime. Env newer than yml => you
# installed/updated something and haven't re-exported.
#
# Advisory only — never blocks, never commits. Exits 0 always.
#   check_env_drift.sh            # human-readable (also prints an all-clear line)
#   check_env_drift.sh --quiet    # print only drifted envs (for the digest)
set -euo pipefail

ENV_DIR="${ENV_DIR:-$HOME/work/scripts/env}"
CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
QUIET=0; [[ "${1:-}" == "--quiet" ]] && QUIET=1

[[ -d "$ENV_DIR" ]] || { echo "check_env_drift: no env dir at $ENV_DIR" >&2; exit 0; }

drift=0
for yml in "$ENV_DIR"/*.yml; do
    [[ -e "$yml" ]] || continue
    env="$(basename "$yml" .yml)"
    meta="$CONDA_ENVS/$env/conda-meta"
    [[ -d "$meta" ]] || continue          # env not built on this node -> can't check, skip
    newest="$(ls -t "$meta"/*.json 2>/dev/null | head -1 || true)"   # || true: head closing the pipe SIGPIPEs ls under pipefail
    [[ -n "$newest" ]] || continue
    if [[ "$newest" -nt "$yml" ]]; then
        echo "[env-drift] '$env' changed since last re-export -> conda env export -n $env > $ENV_DIR/$env.yml"
        drift=$((drift + 1))
    fi
done

if (( QUIET == 0 && drift == 0 )); then
    echo "check_env_drift: all built envs match their committed yml"
fi
exit 0
