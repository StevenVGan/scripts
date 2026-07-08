#!/usr/bin/env bash
# gitall.sh — cross-repo drift detector for the ~/work tree.
# Surfaces un-pushed / dirty repos across nodes (the drift detector referenced
# in PROVISION): walks $WORK, and for each git repo prints one line with its
# working-tree state (clean/DIRTY) and ahead/behind counts vs its upstream.
# Read-only — never touches the index, working tree, or remotes.
#   gitall.sh                # WORK defaults to $HOME/work
#   WORK=<dir> gitall.sh     # scan a different tree
set -euo pipefail

WORK="${WORK:-$HOME/work}"

while IFS= read -r gitdir; do
  r="${gitdir%/.git}"
  if [[ -n "$(git -C "$r" status --porcelain 2>/dev/null)" ]]; then
    state="DIRTY"
  else
    state="clean"
  fi
  if git -C "$r" rev-parse --abbrev-ref '@{u}' >/dev/null 2>&1; then
    ahead="$(git -C "$r" rev-list --count '@{u}..HEAD')"
    behind="$(git -C "$r" rev-list --count 'HEAD..@{u}')"
    sync="ahead ${ahead} behind ${behind}"
  else
    sync="no-upstream"
  fi
  printf '%-6s  %-11s  %s\n' "$state" "$sync" "$r"
done < <(find "$WORK" -maxdepth 3 -type d -name .git -prune)
