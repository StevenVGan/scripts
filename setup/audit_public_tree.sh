#!/usr/bin/env bash
# Refuse to publish institutional/gmail e-mail addresses, UUIDs or secret literals from this PUBLIC repo.
# Location: setup/audit_public_tree.sh
# Usage:  setup/audit_public_tree.sh            # lines a push would add (origin/master..HEAD)   ← run before every push
#         setup/audit_public_tree.sh --staged   # lines staged for the next commit (pre-commit hook mode)
#         setup/audit_public_tree.sh --tree     # every tracked line (baseline sweep)
# Hook (optional, per clone; hooks are not cloned):  ln -s ../../setup/audit_public_tree.sh .git/hooks/pre-commit
# Exit 1 with the offending lines on any hit. Placeholders pass: <COLLECTION_UUID>, you@example.edu, your_password, $VAR.
set -euo pipefail
unset GIT_DIR GIT_WORK_TREE   # git exports GIT_DIR=.git to hooks; resolve the repo from this file's real location instead
cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.."
mode="${1:---push}"
[[ "$0" == *hooks/pre-commit* ]] && mode="--staged"

UUID='[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}'
EMAIL='@(ucsd\.edu|health\.ucsd\.edu|gmail\.com)'
# a secret-looking assignment (any case; the keyword starts a word or follows _) with a literal value; package pins such as tiktoken=0.12 do not match
SECRET_ASSIGN="(^|[^A-Za-z0-9])(secret|token|password|passwd|passphrase|apikey|api_key)=['\"]?[A-Za-z0-9]"
SECRET_PLACEHOLDER="=['\"]?(your_|YOUR_|<|\\\$)"                                     # …unless it is a placeholder
REFRESH='refresh[_]token'

added_lines() {   # "file:line-text" for + lines of a diff
  awk '/^\+\+\+ b\//{f=substr($0,7); next} /^\+/ && !/^\+\+\+/{print f ": " substr($0,2)}'
}
case "$mode" in
  --staged) lines="$(git diff --cached -U0 | added_lines)" ;;
  --tree)   lines="$(git ls-files -z | xargs -0 grep -nHI -E '' -- 2>/dev/null | sed 's/^\([^:]*:[0-9]*\):/\1: /' || true)" ;;
  *)        git rev-parse --verify -q origin/master >/dev/null || { echo "audit: no origin/master; nothing to compare" >&2; exit 0; }
            lines="$( { git diff -U0 origin/master..HEAD | added_lines; git log --format='commit-message(%h): %B' origin/master..HEAD; } )" ;;
esac
n=$( [[ -n "$lines" ]] && grep -c '' <<<"$lines" || echo 0 )
hits="$( {
  printf '%s\n' "$lines" | grep -E "$UUID"      | sed 's/^/UUID     /'
  printf '%s\n' "$lines" | grep -E "$EMAIL"     | sed 's/^/EMAIL    /'
  printf '%s\n' "$lines" | grep -Ei "$SECRET_ASSIGN" | grep -Ev "$SECRET_PLACEHOLDER" | sed 's/^/SECRET   /'
  printf '%s\n' "$lines" | grep -Ei "$REFRESH"  | sed 's/^/TOKEN    /'
} 2>/dev/null | grep -Ev "setup/audit_public_tree.sh:[0-9]+: *(UUID|EMAIL|SECRET_ASSIGN|SECRET_PLACEHOLDER|REFRESH)=" || true )"
if [[ -n "$hits" ]]; then
  echo "audit FAILED ($mode): identifiers or secret literals about to be published:" >&2
  printf '%s\n' "$hits" >&2
  exit 1
fi
echo "audit OK ($mode): $n lines checked, nothing personal or secret found"
