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
# Also checks pip: the mtime heuristic above is structurally blind to pip
# installs (pip writes no conda-meta entry), which is how a broken pyarrow wheel
# hid in `bio` for three months and segfaulted unrelated code through ld.so.
# See scripts/env/README.md "pip installs on this box".
#
# Advisory only — never blocks, never commits. Exits 0 always.
#   check_env_drift.sh            # human-readable (also prints an all-clear line)
#   check_env_drift.sh --quiet    # print only drifted envs (for the digest)
#   PIP_CHECK_ALL=1 check_env_drift.sh   # also sweep envs with no committed yml
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


# ---------------------------------------------------------------------------
# pip-wheel checks (see header). Two signals:
#   (a) pip-BROKEN — a pip package shipping compiled extensions that cannot
#       import. On glibc 2.23 this is the dangerous class: the failed dlopen
#       leaves the linker holding a half-loaded dependency and an unrelated
#       call later segfaults in ld-2.23.so. Always reported.
#   (b) pip-drift — a pip package absent from the yml's pip: section. Only
#       checked when that section came from a real `conda env export` (every
#       entry version-pinned, e.g. bio.yml). Hand-written specs like sc.yml
#       list only DIRECT pip deps while site-packages also holds their
#       transitive deps, so comparing against those is pure noise.
pipbad=0

pip_check_env() {
    local env="$1" yml="${2:-}"
    local envroot="$CONDA_ENVS/$env"
    local py="$envroot/bin/python"
    [[ -x "$py" ]] || return 0
    local sp
    sp="$(ls -d "$envroot"/lib/python*/site-packages 2>/dev/null | head -1 || true)"
    [[ -n "$sp" && -d "$sp" ]] || return 0

    # Parse the yml's pip: section. Entries are indented deeper than the
    # "- pip:" line; blank lines and # comments inside the section are skipped
    # (hand-written specs use them) rather than ending it.
    local declared="" pinned="" has_pip_section=0
    if [[ -n "$yml" ]] && grep -qE '^[[:space:]]*-[[:space:]]*pip:[[:space:]]*$' "$yml"; then
        declared="$(awk '
            /^[[:space:]]*-[[:space:]]*pip:[[:space:]]*$/ { f=1; pi=index($0,"-"); next }
            f && /^[[:space:]]*$/                          { next }
            f && /^[[:space:]]*#/                          { next }
            f && /^[[:space:]]*-[[:space:]]/ {
                     if (index($0,"-") <= pi) exit
                     sub(/^[[:space:]]*-[[:space:]]*/,""); sub(/[[:space:]]*#.*$/,"");
                     print; next }
            f { exit }' "$yml")"
        # Only compare against a section that came from a real `conda env export`
        # (every entry version-pinned). Hand-written specs list only DIRECT pip
        # deps while site-packages also holds their transitive deps, so comparing
        # against those would emit dozens of false positives.
        if [[ -n "$declared" ]]; then
            pinned="$(grep -cv '==' <<<"$declared" || true)"
            [[ "$pinned" == "0" ]] && has_pip_section=1
        fi
        declared="$(sed -E 's/[=<>!~;[].*$//; s/_/-/g' <<<"$declared" | tr '[:upper:]' '[:lower:]')"
    fi

    local dist base name mod norm err
    for dist in "$sp"/*.dist-info; do
        [[ -d "$dist" && -f "$dist/INSTALLER" ]] || continue
        grep -qx 'pip' "$dist/INSTALLER" 2>/dev/null || continue
        base="$(basename "$dist" .dist-info)"
        name="${base%-*}"
        mod=""
        if [[ -f "$dist/top_level.txt" ]]; then
            mod="$(head -1 "$dist/top_level.txt" 2>/dev/null || true)"
        fi
        [[ -n "$mod" ]] || mod="${name//-/_}"
        norm="$(printf '%s' "${name//_/-}" | tr '[:upper:]' '[:lower:]')"

        if (( has_pip_section )) && ! grep -qxF "$norm" <<<"$declared"; then
            echo "[pip-drift]  '$env': pip package '$name' installed but absent from $(basename "$yml") -> conda env export -n $env > $ENV_DIR/$env.yml"
            pipbad=$((pipbad + 1))
        fi

        # only packages with compiled extensions can fail this way
        if [[ -d "$sp/$mod" ]] && find "$sp/$mod" -maxdepth 2 -name '*.so' -print -quit 2>/dev/null | grep -q .; then
            if ! err="$("$py" -c "import $mod" 2>&1)"; then
                echo "[pip-BROKEN] '$env': compiled pip package '$name' cannot import:"
                printf '%s\n' "             $(printf '%s' "$err" | tail -1 | cut -c1-140)"
                echo "             A broken wheel can segfault unrelated code via ld.so. Remove it with:"
                echo "             $py -m pip uninstall -y $name"
                pipbad=$((pipbad + 1))
            fi
        fi
    done
}

if [[ "${PIP_CHECK_ALL:-0}" == "1" ]]; then
    for envroot in "$CONDA_ENVS"/*/; do
        [[ -d "$envroot" ]] || continue
        env="$(basename "$envroot")"
        yml="$ENV_DIR/$env.yml"
        [[ -e "$yml" ]] || yml=""
        pip_check_env "$env" "$yml"
    done
else
    for yml in "$ENV_DIR"/*.yml; do
        [[ -e "$yml" ]] || continue
        env="$(basename "$yml" .yml)"
        [[ -d "$CONDA_ENVS/$env/conda-meta" ]] || continue
        pip_check_env "$env" "$yml"
    done
fi

if (( QUIET == 0 && drift == 0 && pipbad == 0 )); then
    echo "check_env_drift: all built envs match their committed yml; no broken pip wheels"
fi
exit 0
