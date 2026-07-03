#!/usr/bin/env bash
###############################################################################
# fix_project_paths.sh — one-time repair of pre-restructure path rot in the
# 11 bulk-seq project repos (task-6). Guarded per the independent review:
#
#   * EXPLICIT 11-repo allowlist (never a seq/ glob — a glob would also hit
#     4 unreviewed repos: CnR_260212_ZNF432, CnR_260212_ERa_shZNF432,
#     CnR_251204_TOP1cc_HEK, CnR_251204_CTCF_HEK).
#   * Full-line ^BASE= replace to the auto-derive idiom, per-project REL.
#     csRNA is special-cased (../..) because its config lives at
#     scripts/csRNA/0_config.sh (one level deeper than the standard script/).
#   * The 4 INLINE link_fastq.sh wrappers (200207, 240605, 260115, 251215)
#     never define BASE, so a bare ${BASE} swap would abort under set -u —
#     they get a PROJ=dirname/.. derivation + DEST_DIR/MAP_FILE repointed.
#     (The already-correct thin/sourcing wrappers are left untouched.)
#   * Guards: assert exactly one ^BASE= per file; assert the derived BASE
#     resolves to the project root; bash -n after every edit; skip dirty
#     repos; require the digan-mgr-lab remote.
#
# DRY-RUN by default. Pass --apply to write. NEVER commits or pushes —
# per-repo commit+push with sign-off is a separate, manual step.
###############################################################################
set -euo pipefail
WORK="${WORK:-$HOME/work}"
APPLY=0; [[ "${1:-}" == "--apply" ]] && APPLY=1

# entry = seq-relpath | config-relpath | base-REL (dirname/<REL> must == project root)
PROJECTS=(
  "chip/ChIP_GSE59530_ER_p65_MCF7|script/0_config.sh|.."
  "cnr/CnR_200207_ERa_MCF7_Steven|script/0_config.sh|.."
  "cnr/CnR_221220_230207_HA_HepG2-ER3XHA_Steven|script/0_config.sh|.."
  "cnr/CnR_240605_p65_MCF7_TNFa_Likun|script/0_config.sh|.."
  "cnr/CnR_240926_ERa_MCF7_Ctrl_Steven|script/0_config.sh|.."
  "cnr/CnR_260115_ERa_MCF7_OGG1-KD_Priyanka|script/0_config.sh|.."
  "cnr/CnR_260401_ERa_MCF7_OGG1-inhi_Priyanka|script/0_config.sh|.."
  "cnt/CnT_251215_H4K20me1_H3K27me3_MCF7_Steven|script/0_config.sh|.."
  "csrna/csRNA_240216_240701_HepG2-ER3XHA_Steven|scripts/csRNA/0_config.sh|../.."
  "pro/PRO_151113_MCF7_Steven|script/0_config.sh|.."
  "pro/PRO_230222_230321_HepG2-ER3XHA_MCF7_Steven|script/0_config.sh|.."
)
# inline link_fastq.sh wrappers (script/link_fastq.sh) needing PROJ + repoint
LINK_FIX=(
  "cnr/CnR_200207_ERa_MCF7_Steven"
  "cnr/CnR_240605_p65_MCF7_TNFa_Likun"
  "cnr/CnR_260115_ERa_MCF7_OGG1-KD_Priyanka"
  "cnt/CnT_251215_H4K20me1_H3K27me3_MCF7_Steven"
)

[[ ${#PROJECTS[@]} -eq 11 ]] || { echo "FATAL: allowlist != 11 (${#PROJECTS[@]})"; exit 1; }
echo "=== fix_project_paths.sh   mode: $([[ $APPLY -eq 1 ]] && echo '*** APPLY ***' || echo 'DRY-RUN')"
echo "=== WORK=$WORK   (never commits/pushes)"
echo
fail=0

echo "############## PART A — BASE in 0_config.sh (11 repos) ##############"
for entry in "${PROJECTS[@]}"; do
  IFS='|' read -r rel cfgrel baserel <<< "$entry"
  proj="$WORK/seq/$rel"; cfg="$proj/$cfgrel"
  echo "-- $rel   [config: $cfgrel   REL: $baserel]"
  [[ -d "$proj" ]] || { echo "   SKIP: project dir missing"; fail=1; echo; continue; }
  [[ -f "$cfg"  ]] || { echo "   SKIP: config missing";      fail=1; echo; continue; }
  remote=$(git -C "$proj" config --get remote.origin.url 2>/dev/null || echo "")
  [[ "$remote" == *digan-mgr-lab* ]] || { echo "   SKIP: unexpected remote ($remote)"; fail=1; echo; continue; }
  dirty=$(git -C "$proj" status --porcelain 2>/dev/null | wc -l | tr -d ' ')
  [[ "$dirty" == "0" ]] || echo "   WARN: working tree not clean ($dirty change[s]) — inspect before --apply"
  n=$(grep -cE '^[[:space:]]*BASE=' "$cfg" || true)
  [[ "$n" == "1" ]] || { echo "   SKIP: expected 1 BASE= line, found $n"; fail=1; echo; continue; }
  cfgdir="$(cd "$(dirname "$cfg")" && pwd)"; derived="$(cd "$cfgdir/$baserel" && pwd)"; projabs="$(cd "$proj" && pwd)"
  if [[ "$derived" != "$projabs" ]]; then
    echo "   FAIL: derived BASE '$derived' != project root '$projabs' (wrong REL)"; fail=1; echo; continue
  fi
  [[ -e "$proj/samples.tsv" || -d "$proj/data" || -e "$proj/README.md" ]] || echo "   WARN: no samples.tsv/data/README at root"
  old=$(grep -nE '^[[:space:]]*BASE=' "$cfg" | head -1); old="${old#*:}"
  new="BASE=\"\${BASE:-\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")/$baserel\" && pwd)}\""
  echo "   resolves -> $derived  [OK]"
  echo "   OLD: $old"
  echo "   NEW: $new"
  if [[ $APPLY -eq 1 ]]; then
    python3 - "$cfg" "$baserel" <<'PY'
import sys, re
p, rel = sys.argv[1], sys.argv[2]
new = 'BASE="${BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/%s" && pwd)}"\n' % rel
L = open(p).read().splitlines(keepends=True)
idx = [i for i, l in enumerate(L) if re.match(r'^\s*BASE=', l)]
assert len(idx) == 1, "expected exactly one BASE= line, got %d" % len(idx)
L[idx[0]] = new
open(p, 'w').writelines(L)
PY
    bash -n "$cfg" && echo "   APPLIED + bash -n OK"
  fi
  echo
done

echo "############## PART B — inline link_fastq.sh wrappers (4 repos) ##############"
for rel in "${LINK_FIX[@]}"; do
  w="$WORK/seq/$rel/script/link_fastq.sh"
  echo "-- $rel/script/link_fastq.sh"
  [[ -f "$w" ]] || { echo "   SKIP: missing"; fail=1; echo; continue; }
  if grep -qE '^[[:space:]]*(BASE|PROJ)=' "$w"; then
    echo "   SKIP: already defines BASE/PROJ (not an untouched inline wrapper)"; echo; continue
  fi
  grep -qE '^[[:space:]]*DEST_DIR=' "$w" || { echo "   SKIP: no DEST_DIR= line"; fail=1; echo; continue; }
  echo "   OLD:"; grep -nE '^[[:space:]]*(DEST_DIR|MAP_FILE)=' "$w" | sed 's/^/     /'
  echo "   NEW (RAW_DIR left unchanged):"
  echo "     PROJ=\"\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")/..\" && pwd)\""
  echo "     DEST_DIR=\"\${DEST_DIR:-\$PROJ/data}\""
  echo "     MAP_FILE=\"\${MAP_FILE:-\$PROJ/link_sample.tsv}\""
  if [[ $APPLY -eq 1 ]]; then
    python3 - "$w" <<'PY'
import sys, re
p = sys.argv[1]
L = open(p).read().splitlines(keepends=True)
out = []; inserted = False
for l in L:
    if re.match(r'^\s*DEST_DIR=', l) and not inserted:
        out.append('PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # project root = parent of script/\n')
        out.append('DEST_DIR="${DEST_DIR:-$PROJ/data}"\n')
        inserted = True
        continue
    if re.match(r'^\s*MAP_FILE=', l):
        out.append('MAP_FILE="${MAP_FILE:-$PROJ/link_sample.tsv}"\n')
        continue
    out.append(l)
assert inserted, "no DEST_DIR= line found to anchor the PROJ insertion"
open(p, 'w').writelines(out)
PY
    bash -n "$w" && echo "   APPLIED + bash -n OK"
  fi
  echo
done

if [[ $fail -eq 0 ]]; then echo "=== ALL CHECKS PASSED"; else echo "=== SOME CHECKS FAILED (see SKIP/FAIL above)"; fi
echo "=== Next (manual, per repo): review diff -> commit -> push, each with sign-off. This script never commits."
[[ $fail -eq 0 ]] || exit 1
