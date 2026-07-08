#!/usr/bin/env bash
# Mirror the work_brain repo (working tree + .git) to a second location as a
# disk-failure safety net for state not yet pushed to GitHub — this captures
# BOTH committed-but-unpushed commits (in .git) AND uncommitted working-tree
# edits. GitHub remains the off-site backup of *pushed* history; this closes the
# gap between pushes against a local disk failure on the plan-home node.
#
# It is a MIRROR, not versioned history (git is the history). Run from cron on
# the plan-home node, e.g. hourly:
#   0 * * * * $HOME/work/scripts/setup/backup_work_brain.sh >> $HOME/.node_guard/work_brain_backup.log 2>&1
set -euo pipefail

SRC="${WORK_BRAIN:-$HOME/work/work_brain}"
DEST="${WORK_BRAIN_BACKUP:-/mnt/share/archive/bkup/work_brain_mirror}"

[[ -d "$SRC" ]] || { echo "backup_work_brain: source $SRC not found" >&2; exit 1; }

# Bail quietly (exit 0) if the target can't be created — NAS unmounted, or the
# path not writable — so a cron run doesn't error-spam. Point WORK_BRAIN_BACKUP
# at a writable location (the /mnt/share default assumes a writable NAS dir).
if ! mkdir -p "$DEST" 2>/dev/null; then
    echo "backup_work_brain: cannot create $DEST (NAS unmounted or not writable) — set WORK_BRAIN_BACKUP to a writable path; skipping" >&2
    exit 0
fi
rsync -a --delete "$SRC"/ "$DEST"/
echo "backup_work_brain: mirrored $SRC -> $DEST at $(date '+%Y-%m-%d %H:%M:%S')"
