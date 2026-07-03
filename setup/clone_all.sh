#!/usr/bin/env bash
# clone_all.sh — clone the shared spine onto a new node (HTTPS only).
# Default: scripts/ (public spine) + work_brain/ (private durable Claude brain).
# Add the standalone tool repos with flags:
#   clone_all.sh [--all] [--tf-atlas] [--pwm] [--idr]
# WORK=<dir> to place the tree somewhere other than $HOME/work.
set -euo pipefail

ALL=0; TF=0; PWM=0; IDR=0
for a in "$@"; do case "$a" in
  --all) ALL=1 ;; --tf-atlas) TF=1 ;; --pwm) PWM=1 ;; --idr) IDR=1 ;;
  -h|--help) grep '^# ' "$0" | sed 's/^# //'; exit 0 ;;
  *) echo "unknown flag: $a (see --help)"; exit 2 ;;
esac; done

WORK="${WORK:-$HOME/work}"; mkdir -p "$WORK"; cd "$WORK"

clone() {  # url dir
  if [[ -d "$2/.git" ]]; then echo "skip (exists): $2"; else echo "clone: $2"; git clone "$1" "$2"; fi
}

clone https://github.com/StevenVGan/scripts.git scripts
# PRIVATE durable Claude brain — non-fatal so scripts/ still lands if it's not created/pushed yet
clone https://github.com/StevenVGan/work_brain.git work_brain \
  || echo "warn: work_brain not cloned — create+push the PRIVATE repo first (see bootstrap.md)"
if [[ $ALL -eq 1 || $TF  -eq 1 ]]; then clone https://github.com/StevenVGan/tf_atlas.git tf_atlas; fi
if [[ $ALL -eq 1 || $PWM -eq 1 ]]; then
  clone https://github.com/StevenVGan/PWM-motif-analysis.git "PWM motif analysis"
  clone https://github.com/StevenVGan/PWM_motif_analysis_v2.git PWM_motif_analysis_v2
fi
if [[ $ALL -eq 1 || $IDR -eq 1 ]]; then clone https://github.com/digan-mgr-lab/IDR_analysis.git IDR_analysis; fi

echo
echo "Next: rebuild envs  ->  mamba env create -f scripts/env/bio.yml ; mamba env create -f scripts/env/sc.yml"
echo "Then scaffold a project  ->  scripts/setup/new_project.sh <assay> <name>"
echo "See scripts/setup/bootstrap.md for the full new-node runbook."
