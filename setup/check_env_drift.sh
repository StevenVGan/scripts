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
#   (b) pip-drift — a pip package installed but absent from the yml. Only
#       checked for EXPORT-style ymls (`conda env export` output, recognised by
#       a `prefix:` line OR any `name=version=build` pin -- pwm2.yml is an
#       export whose `prefix:` was stripped), which claim to enumerate
#       everything. Hand-written specs (sc.yml, meth.yml) list only DIRECT pip
#       deps while site-packages also holds their transitive deps, so comparing
#       against those is noise.
pipbad=0

pip_check_env() {
    local env="$1" yml="${2:-}"
    local envroot="$CONDA_ENVS/$env"
    local py="$envroot/bin/python"
    [[ -x "$py" ]] || return 0
    local sp
    sp="$(ls -d "$envroot"/lib/python*/site-packages 2>/dev/null | head -1 || true)"
    [[ -n "$sp" && -d "$sp" ]] || return 0

    # Is this yml an EXPORT (`conda env export`) or a hand-written spec?
    # Only an export claims to enumerate everything, so only there is "installed
    # but undeclared" real drift. Markers: a trailing `prefix:` line, OR any
    # dependency pinned `name=version=build` -- exports pin every one of them
    # (bio 327, rna 168, pwm2 136), specs pin none (sc.yml, meth.yml both 0).
    # Both tests are needed: pwm2.yml is an export whose `prefix:` line was
    # stripped, so a prefix-only test would misfile it as a spec and skip it.
    # Key off this, NOT off a pip: section existing -- an export written when the
    # env had no pip packages has no pip: section at all, and every pip package
    # installed since is exactly the drift to catch (this is how `primer` hid
    # two undeclared pip packages).
    local declared="" is_export=0
    if [[ -n "$yml" ]] && { grep -qE '^prefix:' "$yml" \
                            || grep -qE '^[[:space:]]*-[[:space:]][^ ]+=[^ ]+=' "$yml"; }; then
        is_export=1
        declared="$(awk '
            /^[[:space:]]*-[[:space:]]*pip:[[:space:]]*$/ { f=1; pi=index($0,"-"); next }
            f && /^[[:space:]]*$/                          { next }
            f && /^[[:space:]]*#/                          { next }
            f && /^[[:space:]]*-[[:space:]]/ {
                     if (index($0,"-") <= pi) exit
                     sub(/^[[:space:]]*-[[:space:]]*/,""); sub(/[[:space:]]*#.*$/,"");
                     print; next }
            f { exit }' "$yml" \
            | sed -E 's/[=<>!~;[].*$//; s/_/-/g' | tr '[:upper:]' '[:lower:]')"
    fi

    local dist base name norm cands err
    for dist in "$sp"/*.dist-info; do
        [[ -d "$dist" && -f "$dist/INSTALLER" ]] || continue
        grep -qx 'pip' "$dist/INSTALLER" 2>/dev/null || continue
        base="$(basename "$dist" .dist-info)"
        name="${base%-*}"
        norm="$(printf '%s' "${name//_/-}" | tr '[:upper:]' '[:lower:]')"

        if (( is_export )) && ! grep -qxF "$norm" <<<"$declared"; then
            echo "[pip-drift]  '$env': pip package '$name' installed but absent from $(basename "$yml") -> conda env export -n $env > $ENV_DIR/$env.yml"
            pipbad=$((pipbad + 1))
        fi

        # Only a dist shipping compiled extensions can fail the way a bad
        # manylinux wheel does. Read RECORD, not the filesystem: a top_level.txt
        # + directory test misses root-level .so files (pillow_heif ships
        # _pillow_heif.cpython-311-*.so directly in site-packages).
        [[ -f "$dist/RECORD" ]] || continue
        grep -qE '\.so(\.[0-9]+)*,' "$dist/RECORD" || continue

        # Candidate import names: EVERY line of top_level.txt (not just the
        # first), plus top-level names inferred from RECORD -- packages with an
        # __init__.py and root-level .py/.so modules. Needed because many dists
        # ship no top_level.txt and their module name differs from the dist name
        # (scikit_image -> skimage, protobuf -> google).
        cands="$( { [[ -f "$dist/top_level.txt" ]] && tr -d "\r" < "$dist/top_level.txt"
                    cut -d, -f1 "$dist/RECORD" | awk -F/ '
                        $0 ~ /\.dist-info\// { next }
                        NF>1 && $2=="__init__.py" { print $1; next }
                        NF==1 { f=$1; sub(/\.(py|so)$/,"",f); sub(/\.cpython-[^.]*$/,"",f); print f }' ; } \
                  | sed 's/[[:space:]]//g' | grep -E '^[A-Za-z_][A-Za-z0-9_]*$' | sort -u || true)"
        [[ -n "$cands" ]] || continue

        # One interpreter start per dist; it tries each candidate and stops at the
        # first that imports. A candidate whose OWN name is simply wrong raises
        # ModuleNotFoundError for itself and is skipped; any other failure -- a
        # missing GLIBC symbol, a missing dependency -- is a real break.
        if ! err="$("$py" -c '
import importlib, sys
for m in sys.argv[1:]:
    try:
        importlib.import_module(m)
    except ModuleNotFoundError as e:
        if getattr(e, "name", None) == m:
            continue
        print("%s: %s" % (type(e).__name__, e)); sys.exit(1)
    except BaseException as e:
        print("%s: %s" % (type(e).__name__, e)); sys.exit(1)
    else:
        sys.exit(0)
sys.exit(0)
' $cands 2>&1)"; then
            echo "[pip-BROKEN] '$env': compiled pip package '$name' cannot import:"
            printf '%s\n' "             $(printf '%s' "$err" | tail -1 | cut -c1-140)"
            echo "             A broken wheel can segfault unrelated code via ld.so. Remove it with:"
            echo "             $py -m pip uninstall -y $name"
            pipbad=$((pipbad + 1))
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
