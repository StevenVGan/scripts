#!/usr/bin/env bash
###############################################################################
# new_project.sh — scaffold a new bulk-seq project per scripts/CONVENTIONS.md
#
# Stamps everything that is pure structure/boilerplate; leaves the biology
# (sample rows, IP/control pairings, raw-FASTQ map, README prose) as clearly
# marked TODOs. Nothing is committed — it `git init`s and stages, then prints
# the review + next steps so you can inspect before the first commit.
#
# Usage:
#   new_project.sh [options] <pipeline> <project_name>
#
#   <pipeline>       cnr | csrna | pro | atac   (which scripts/pipeline/<dir> to copy)
#   <project_name>   e.g. CnR_260501_ERa_MCF7_E2_Steven   (Assay_date_target_cell[_treat]_person)
#
# Options:
#   --dest <subdir>   place under seq/<subdir>/ (default: <pipeline>; use for
#                     chip/cnt reusing the cnr pipeline, e.g. --dest cnt)
#   --label <text>    README "Assay:" label (default: mapped from pipeline)
#   --genome <g>      genome (default: hg38)
#   --se              single-end (default: paired-end)
#   --owner <name>    README owner/submitter
#   --upstream <txt>  README upstream source (GEO accession / IGM submission / collaborator)
#   --remote <url>    git remote to add (HTTPS; not pushed — printed for you to push)
#   --dry-run         print the plan, create nothing
#   -h | --help
#
# Examples:
#   new_project.sh cnr CnR_260501_ERa_MCF7_E2_Steven --owner Steven --upstream "IGM 260501"
#   new_project.sh cnr CnT_260501_H3K27me3_MCF7_Steven --dest cnt --label "CUT&Tag"
#   new_project.sh csrna csRNA_260601_HepG2_Steven --se
###############################################################################
set -euo pipefail

die() { echo "[ERROR] $*" >&2; exit 1; }
info() { echo "[new_project] $*"; }
# escape a string for use as a sed replacement (\, &, and the | delimiter)
esc() { printf '%s' "$1" | sed -e 's/[\\&|]/\\&/g'; }

# ---- locate scripts/ (this file lives in scripts/setup/) --------------------
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMPLATES="${SCRIPTS_DIR}/pipeline/templates"
WORK_ROOT="${WORK_ROOT:-$(cd "${SCRIPTS_DIR}/.." && pwd)}"   # default: parent of scripts/ (= ~/work)
[[ -d "$TEMPLATES" ]] || die "templates dir not found: $TEMPLATES"

# ---- defaults + arg parsing -------------------------------------------------
GENOME="hg38"; SE="0"; DEST=""; LABEL=""; OWNER="TODO"; UPSTREAM="TODO"; REMOTE=""; DRY=0
PIPELINE=""; PROJECT=""
positional=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest)     DEST="$2"; shift 2 ;;
    --label)    LABEL="$2"; shift 2 ;;
    --genome)   GENOME="$2"; shift 2 ;;
    --se)       SE="1"; shift ;;
    --owner)    OWNER="$2"; shift 2 ;;
    --upstream) UPSTREAM="$2"; shift 2 ;;
    --remote)   REMOTE="$2"; shift 2 ;;
    --dry-run)  DRY=1; shift ;;
    -h|--help)  sed -n '2,32p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    -*)         die "unknown option: $1" ;;
    *)          positional+=("$1"); shift ;;
  esac
done
[[ ${#positional[@]} -eq 2 ]] || die "expected <pipeline> <project_name>; got ${#positional[@]} positional args. See --help."
PIPELINE="${positional[0]}"; PROJECT="${positional[1]}"

case "$PIPELINE" in
  cnr|csrna|pro|atac) : ;;
  *) die "pipeline must be one of: cnr csrna pro atac (got '$PIPELINE')" ;;
esac
PIPELINE_SRC="${SCRIPTS_DIR}/pipeline/${PIPELINE}"
[[ -d "$PIPELINE_SRC" ]] || die "pipeline source not found: $PIPELINE_SRC"

DEST="${DEST:-$PIPELINE}"
LAYOUT=$([[ "$SE" == "1" ]] && echo SE || echo PE)
if [[ -z "$LABEL" ]]; then
  case "$PIPELINE" in
    cnr)   LABEL="CUT&RUN" ;;
    csrna) LABEL="csRNA" ;;
    pro)   LABEL="PRO-seq" ;;
    atac)  LABEL="ATAC" ;;
  esac
fi
DATE="$(date +%F)"
PROJDIR="${WORK_ROOT}/seq/${DEST}/${PROJECT}"

info "pipeline   : $PIPELINE  (src: $PIPELINE_SRC)"
info "project    : seq/${DEST}/${PROJECT}"
info "genome/SE  : ${GENOME} / ${LAYOUT}"
info "dest path  : $PROJDIR"
[[ $DRY -eq 1 ]] && { info "(dry-run — no changes made)"; exit 0; }
[[ -e "$PROJDIR" ]] && die "destination already exists: $PROJDIR (refusing to overwrite)"

# ---- 1. directory tree ------------------------------------------------------
info "creating directory tree"
mkdir -p "$PROJDIR"/{data,cleandata,align/{bam,track,tags},peaks/{macs3,homer},multiqc,logs,scratch,script/local,figures,qc,analysis}

cd "$PROJDIR"

# ---- 2. pipeline step scripts (FLAT into script/; untracked by gitignore) ---
info "copying pipeline scripts (flat)"
cp "$PIPELINE_SRC"/*.sh script/
chmod +x script/*.sh

# bake GENOME / SE into the copied 0_config.sh (BASE auto-derives — never touched)
if [[ "$GENOME" != "hg38" ]]; then
  sed -i -E "s|^GENOME=.*|GENOME=\"\${GENOME:-${GENOME}}\"|" script/0_config.sh || true
fi
if [[ "$SE" == "1" ]]; then
  sed -i -E "s|^SE=.*|SE=\"\${SE:-1}\"|" script/0_config.sh || true
fi

# ---- 3. link_fastq.sh wrapper (exec's shared prep — no inlined logic) -------
info "writing link_fastq.sh wrapper"
cat > script/link_fastq.sh <<'WRAP'
#!/usr/bin/env bash
# Project link_fastq wrapper — exports paths, exec's the shared prep script.
# Logic lives in scripts/pipeline/tools/prep/link_fastq.sh — do NOT inline it here.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${BASE:-$(cd "$HERE/.." && pwd)}"
PREP="${PREP:-$HOME/work/scripts/pipeline/tools/prep/link_fastq.sh}"

export RAW_DIR="${RAW_DIR:-/PATH/TO/raw_seq/RUN}"   # TODO: raw FASTQ source dir
export DEST_DIR="${DEST_DIR:-$BASE/data}"           # derived from BASE — do not point at a test dir
export MAP_FILE="${MAP_FILE:-$BASE/link_sample.tsv}"

exec "$PREP"
WRAP
chmod +x script/link_fastq.sh

# ---- 4. .gitignore + metadata templates -------------------------------------
info "emitting .gitignore + metadata templates"
cp "$TEMPLATES/gitignore"            .gitignore
cp "$TEMPLATES/samples.tsv"          samples.tsv
cp "$TEMPLATES/peakcall_groups.tsv"  peakcall_groups.tsv
cp "$TEMPLATES/link_sample.tsv"      link_sample.tsv
if [[ "$PIPELINE" == "csrna" ]]; then
  cp "$TEMPLATES/tss_groups.tsv"     tss_groups.tsv    # csRNA only; PRO's 4.3 uses no groups file
fi

printf '# Milestone figures\n\n| File | Generated by | Date | Notes |\n|---|---|---|---|\n' > figures/manifest.md

# ---- 5. README from template (token substitution) ---------------------------
info "writing README.md"
sed -e "s|__PROJECT__|$(esc "$PROJECT")|g" \
    -e "s|__ASSAY__|$(esc "$LABEL")|g" \
    -e "s|__ASSAY_DIR__|$(esc "$PIPELINE")|g" \
    -e "s|__GENOME__|$(esc "$GENOME")|g" \
    -e "s|__LAYOUT__|$(esc "$LAYOUT")|g" \
    -e "s|__OWNER__|$(esc "$OWNER")|g" \
    -e "s|__DATE__|$(esc "$DATE")|g" \
    -e "s|__UPSTREAM__|$(esc "$UPSTREAM")|g" \
    "$TEMPLATES/README.project.md" > README.md

# ---- 6. git init + stage (NOT committed) ------------------------------------
info "git init + stage (not committed)"
git init -q
git symbolic-ref HEAD refs/heads/main   # unborn 'main' (git 2.7.4 has no init -b)
[[ -n "$REMOTE" ]] && git remote add origin "$REMOTE"
git add -A

# ---- done -------------------------------------------------------------------
cat <<DONE

[new_project] scaffolded: $PROJDIR

Stamped (structure/boilerplate — done):
  - dir tree, pipeline scripts (flat, untracked), .gitignore (one-per-line)
  - samples.tsv / peakcall_groups.tsv / link_sample.tsv$([[ "$PIPELINE" == csrna ]] && echo " / tss_groups.tsv") headers
  - README.md (Assay=$LABEL, Genome=$GENOME, $LAYOUT), figures/manifest.md
  - link_fastq.sh wrapper (exec's shared prep; DEST_DIR/MAP_FILE derived from BASE)
  - 0_config.sh: BASE auto-derives$([[ "$GENOME" != hg38 ]] && echo "; GENOME baked")$([[ "$SE" == 1 ]] && echo "; SE=1 baked")

Now supply the biology (TODO — the parts a scaffolder can't):
  1. link_sample.tsv  — raw_prefix -> sample_name map; set RAW_DIR in script/link_fastq.sh
  2. samples.tsv      — one row per sample (target/condition/replicate/fastq paths)
  3. peakcall_groups.tsv — IP/control pairings + TF|Histone
  4. README.md        — fill "What this project is" + Owner/Upstream if left TODO
  5. run:  cd script && ./link_fastq.sh && ./run_all.sh   # references.tsv auto-emitted by 5_qc.sh

Then review + commit:
  git -C "$PROJDIR" status
  git -C "$PROJDIR" commit -m "scaffold ${PROJECT}"
$([[ -n "$REMOTE" ]] && echo "  git -C \"$PROJDIR\" push -u origin main   # remote: $REMOTE" || echo "  # add a remote when ready: gh repo create ... ; git remote add origin <https-url>")
DONE
