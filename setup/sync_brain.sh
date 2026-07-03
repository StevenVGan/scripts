#!/usr/bin/env bash
# sync_brain.sh — sync the DURABLE Claude memory between this node's live memory
# dir and the git-tracked carry-bundle in the PRIVATE work_brain repo
# ($HOME/work/work_brain/memory/).
#
#   push : live  -> bundle   (feedback_* + user_* + reference_*, minus node-specific)
#   pull : bundle -> live     (durable set; leaves box-local project_* untouched)
#
# Filter (memory `type` taxonomy + PROVISION App-G):
#   feedback_*, user_*, reference_*  -> synced EXCEPT the node-specific facts below
#                                       (facts describing one machine, not the lab)
#   project_*                        -> never collected by push (box-local findings);
#                                       the 2 seeded project notes are pull-only.
set -euo pipefail
# durable-memory bundle lives in the PRIVATE work_brain repo (not public scripts/)
BUNDLE="${BRAIN_REPO:-$HOME/work/work_brain}/memory"
# live Claude memory dir; override with LIVE=... . Default = the per-project dir
# Claude Code uses for $HOME/work (path with '/' -> '-').
LIVE="${LIVE:-$HOME/.claude/projects/$(printf '%s' "$HOME/work" | sed 's#/#-#g')/memory}"

# memories that describe ONE machine (not the lab) — carried once as re-validate
# seeds, never auto-synced. Matches MEMORY.md's re-validate header + PROVISION App-G.
NODE_SPECIFIC="reference_git_https_only reference_gh_binary_location reference_genome_fastas feedback_conda_env_build_linux01 feedback_cpu_cap_16"
is_node_specific() { case " $NODE_SPECIFIC " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

case "${1:-}" in
  push)
    [[ -d "$LIVE" ]] || { echo "live memory dir not found: $LIVE (set LIVE=...)"; exit 1; }
    mkdir -p "$BUNDLE"
    for f in "$LIVE"/feedback_*.md "$LIVE"/user_*.md "$LIVE"/reference_*.md; do
      [[ -e "$f" ]] || continue
      is_node_specific "$(basename "$f" .md)" || cp "$f" "$BUNDLE/"
    done
    echo "pushed durable memories -> $BUNDLE"
    echo "NOTE: regenerate $BUNDLE/MEMORY.md to index the new set, review the diff,"
    echo "      then commit (per-commit sign-off; no Co-Authored-By trailer)."
    ;;
  pull)
    [[ -d "$BUNDLE" ]] || { echo "bundle not found: $BUNDLE"; exit 1; }
    mkdir -p "$LIVE"
    for f in "$BUNDLE"/*.md; do
      [[ "$(basename "$f")" == "MEMORY.md" ]] && continue
      cp "$f" "$LIVE/"
    done
    echo "pulled durable memories -> $LIVE (box-local project_* left untouched)."
    echo "NOTE: merge $BUNDLE/MEMORY.md lines into $LIVE/MEMORY.md by hand so this"
    echo "      node's own project_* index lines are not clobbered."
    ;;
  *) echo "usage: sync_brain.sh {push|pull}   (LIVE=<memory dir> to override)"; exit 2 ;;
esac
