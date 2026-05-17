#!/usr/bin/env bash
# scripts/maintenance/archive_inactive.sh — see README.md for the full process.
#
# Subcommands:
#   inventory               build inventory.tsv (sizes per cleanup project; no deletions)
#   run                     archive ALL cleanup-status projects from projects.tsv
#   run <project_path>      archive a single project (smoke test)
#   audit                   compare current sizes vs inventory.tsv
#
# Per-project actions on cleanup-status rows only:
#   1. BAM -> CRAM (lossless), verify (read-count + samtools quickcheck), then delete the BAM + .bai
#   2. delete align/tags/   (regenerable via HOMER makeTagDirectory)
#   3. delete cleandata/, cleandata2/   (regenerable via 1_trim_qc from raw_seq/)
# Bigwigs (align/track/), peaks/, multiqc/, data/, script/, logs/ are NEVER touched.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HOME/work}"
PROJECTS_TSV="${PROJECTS_TSV:-$HERE/projects.tsv}"
INVENTORY_TSV="${INVENTORY_TSV:-$HERE/inventory.tsv}"
AUDIT_TSV="${AUDIT_TSV:-$HERE/audit.tsv}"
LOG_DIR="${LOG_DIR:-$HERE/logs}"
mkdir -p "$LOG_DIR"

# Reference FASTAs by genome — override via env if paths change
REF_HG38="${REF_HG38:-/mnt/share/archive/bkup/ref/genome/hg38/hg38.fa}"
REF_HG19="${REF_HG19:-/mnt/share/archive/bkup/ref/genome/hg19/hg19.fa}"
REF_MM10="${REF_MM10:-/mnt/share/archive/bkup/ref/genome/mm10/mm10.fa}"
DEFAULT_GENOME="${DEFAULT_GENOME:-}"   # used only if auto-detect fails; set to e.g. hg38 to assume

THREADS="${THREADS:-8}"

log()  { printf '[%(%FT%T)T] %s\n' -1 "$*" >&2; }
fail() { log "ERROR: $*"; exit 1; }
du_bytes() { du -sb "$1" 2>/dev/null | awk '{print $1; ok=1} END {if (!ok) print 0}'; }

ref_for() {
  case "${1:-}" in
    hg38) echo "$REF_HG38" ;;
    hg19) echo "$REF_HG19" ;;
    mm10) echo "$REF_MM10" ;;
    *)    echo "" ;;
  esac
}

find_bams() {
  # Standard layout: align/bam/*.bam ; older projects: align/*.bam
  local proj="$1"
  { find "$proj/align/bam" -maxdepth 1 -type f -name "*.bam" 2>/dev/null
    find "$proj/align"     -maxdepth 1 -type f -name "*.bam" 2>/dev/null
  } | sort -u
}

detect_genome() {
  local proj="$1" cfg="$proj/script/0_config.sh" g=""
  if [[ -f "$cfg" ]]; then
    g=$(awk -F= '/^GENOME=/{gsub(/[" \047]/,"",$2); print $2; exit}' "$cfg")
  fi
  if [[ -z "$g" ]] && command -v samtools >/dev/null 2>&1; then
    local bam len
    bam=$(find_bams "$proj" | head -1)
    if [[ -n "$bam" ]]; then
      len=$(samtools view -H "$bam" 2>/dev/null \
        | awk -F'[\t:]' '/^@SQ/ && $3=="chr1"{print $5; exit}')
      case "$len" in
        248956422) g=hg38 ;;
        249250621) g=hg19 ;;
        195471971) g=mm10 ;;
      esac
    fi
  fi
  [[ -z "$g" ]] && g="$DEFAULT_GENOME"
  echo "$g"
}

# ---------------- inventory ----------------
cmd_inventory() {
  printf 'project\tstatus\tgenome\tbam_n\tbam_bytes\ttags_bytes\tcleandata_bytes\ttrack_bytes\tpeaks_bytes\ttotal_bytes\n' > "$INVENTORY_TSV"
  while IFS=$'\t' read -r status path notes; do
    [[ "$status" == "cleanup" ]] || continue
    local proj="$WORK/$path"
    if [[ ! -d "$proj" ]]; then log "missing: $path"; continue; fi

    local genome bam_n bam_b tags_b c1 c2 clean_b track_b peaks_b total_b
    genome=$(detect_genome "$proj"); genome=${genome:-unknown}

    bam_n=$(find_bams "$proj" | grep -c . || true)
    bam_b=0
    while IFS= read -r f; do
      [[ -z "$f" ]] && continue
      bam_b=$(( bam_b + $(stat -c%s "$f" 2>/dev/null || echo 0) ))
    done < <(find_bams "$proj")

    tags_b=$(du_bytes  "$proj/align/tags");  tags_b=${tags_b:-0}
    c1=$(du_bytes      "$proj/cleandata");   c1=${c1:-0}
    c2=$(du_bytes      "$proj/cleandata2");  c2=${c2:-0}
    clean_b=$((c1 + c2))
    track_b=$(du_bytes "$proj/align/track"); track_b=${track_b:-0}
    peaks_b=$(du_bytes "$proj/peaks");       peaks_b=${peaks_b:-0}
    total_b=$(du_bytes "$proj");             total_b=${total_b:-0}

    printf '%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' \
      "$path" "$status" "$genome" \
      "$bam_n" "$bam_b" "$tags_b" "$clean_b" "$track_b" "$peaks_b" "$total_b" \
      >> "$INVENTORY_TSV"
  done < <(tail -n +2 "$PROJECTS_TSV")

  log "wrote $INVENTORY_TSV"
  awk -F'\t' 'NR>1 {bam+=$5; tags+=$6; clean+=$7; tot+=$10}
              END {printf "Cleanup-set baseline: %.1f GB total; bam=%.1f GB tags=%.1f GB cleandata=%.1f GB\n",
                          tot/2^30, bam/2^30, tags/2^30, clean/2^30}' "$INVENTORY_TSV"
}

# ---------------- archive a single project ----------------
archive_one_body() {
  local proj="$1"
  log "==> $proj"

  local genome
  genome=$(detect_genome "$proj")
  log "genome detected: ${genome:-NONE}"

  # 1) BAM -> CRAM
  local n_bams
  n_bams=$(find_bams "$proj" | grep -c . || true)
  if [[ "$n_bams" -gt 0 ]]; then
    local ref
    ref=$(ref_for "$genome")
    if [[ -z "$ref" || ! -f "$ref" ]]; then
      log "  WARN: no reference for genome='${genome:-?}' (looked for: ${ref:-<none>}). Skipping BAM->CRAM."
    elif ! command -v samtools >/dev/null 2>&1; then
      log "  WARN: samtools not in PATH; activate the bio env. Skipping BAM->CRAM."
    else
      while IFS= read -r bam; do
        [[ -z "$bam" ]] && continue
        local cram="${bam%.bam}.cram"
        if [[ -f "$cram" ]]; then
          log "  cram exists, skipping bam: $(basename "$bam")"
          continue
        fi
        log "  cram: $(basename "$bam")"
        if ! samtools view -@ "$THREADS" -C -T "$ref" -o "$cram" "$bam"; then
          log "  FAIL view; keeping bam"; rm -f "$cram"; continue
        fi
        if ! samtools index -@ "$THREADS" "$cram"; then
          log "  FAIL index; keeping bam"; rm -f "$cram" "$cram.crai"; continue
        fi
        local nb nc
        nb=$(samtools view -@ "$THREADS" -c "$bam" 2>/dev/null)
        nc=$(samtools view -@ "$THREADS" -c -T "$ref" "$cram" 2>/dev/null)
        if [[ "$nb" != "$nc" || -z "$nb" || -z "$nc" ]]; then
          log "  VERIFY FAIL: bam=$nb cram=$nc; keeping bam"
          rm -f "$cram" "$cram.crai"; continue
        fi
        if ! samtools quickcheck -v "$cram" 2>/dev/null; then
          log "  quickcheck FAIL; keeping bam"
          rm -f "$cram" "$cram.crai"; continue
        fi
        rm -f "$bam" "${bam}.bai" "${bam%.bam}.bai"
        log "  ok ($nb reads, $(numfmt --to=iec --suffix=B "$(stat -c%s "$cram")"))"
      done < <(find_bams "$proj")
    fi
  else
    log "  no BAMs to compress"
  fi

  # 2) tag dirs
  if [[ -d "$proj/align/tags" ]]; then
    log "  removing align/tags ($(du -sh "$proj/align/tags" | awk '{print $1}'))"
    rm -rf "$proj/align/tags"
  fi

  # 3) trimmed fastqs
  for d in cleandata cleandata2; do
    if [[ -d "$proj/$d" ]]; then
      log "  removing $d ($(du -sh "$proj/$d" | awk '{print $1}'))"
      rm -rf "$proj/$d"
    fi
  done

  log "<== done"
}

archive_one() {
  local proj="$1"
  [[ -d "$proj" ]] || fail "no such project: $proj"
  local plog="$LOG_DIR/$(basename "$proj").$(date +%Y%m%d-%H%M%S).log"
  ( archive_one_body "$proj" ) 2>&1 | tee -a "$plog"
}

# ---------------- run (batch or single) ----------------
cmd_run() {
  if [[ $# -gt 0 ]]; then
    local arg="$1" proj="$1"
    [[ -d "$proj" ]] || proj="$WORK/$arg"
    [[ -d "$proj" ]] || fail "no such project: $arg"
    archive_one "$proj"
    return
  fi
  local n_done=0 n_fail=0
  while IFS=$'\t' read -r status path notes; do
    [[ "$status" == "cleanup" ]] || continue
    local proj="$WORK/$path"
    if [[ ! -d "$proj" ]]; then log "missing: $path (skip)"; continue; fi
    if archive_one "$proj"; then n_done=$((n_done+1)); else n_fail=$((n_fail+1)); log "FAILED: $path"; fi
  done < <(tail -n +2 "$PROJECTS_TSV")
  log "batch finished: $n_done ok, $n_fail failed"
}

# ---------------- audit ----------------
cmd_audit() {
  [[ -f "$INVENTORY_TSV" ]] || fail "no inventory at $INVENTORY_TSV — run 'inventory' first"
  printf 'project\tbefore_bytes\tafter_bytes\tsaved_bytes\tsaved_pct\n' > "$AUDIT_TSV"
  while IFS=$'\t' read -r project status genome bam_n bam_b tags_b clean_b track_b peaks_b total_b; do
    [[ "$project" == "project" ]] && continue
    local proj="$WORK/$project"
    [[ -d "$proj" ]] || continue
    local after saved pct=0
    after=$(du_bytes "$proj"); after=${after:-0}
    saved=$((total_b - after))
    [[ "$total_b" -gt 0 ]] && pct=$(( saved*100/total_b ))
    printf '%s\t%d\t%d\t%d\t%d\n' "$project" "$total_b" "$after" "$saved" "$pct" >> "$AUDIT_TSV"
  done < "$INVENTORY_TSV"
  log "wrote $AUDIT_TSV"
  awk -F'\t' 'NR>1 {b+=$2; a+=$3; s+=$4}
              END {printf "Total before: %.1f GB\nTotal after:  %.1f GB\nSaved:        %.1f GB (%.1f%%)\n",
                          b/2^30, a/2^30, s/2^30, (b? s*100/b: 0)}' "$AUDIT_TSV"
}

# ---------------- dispatch ----------------
case "${1:-}" in
  inventory) shift; cmd_inventory "$@" ;;
  run)       shift; cmd_run "$@" ;;
  audit)     shift; cmd_audit "$@" ;;
  *) cat <<EOF >&2
usage: $(basename "$0") {inventory|run|audit} [project_path]

  inventory               build inventory.tsv (sizes per cleanup project; no deletions)
  run                     archive ALL cleanup-status projects from projects.tsv
  run <project_path>      archive a single project (smoke test)
  audit                   compare current sizes vs inventory.tsv

env overrides: REF_HG38, REF_HG19, REF_MM10, DEFAULT_GENOME, THREADS,
               PROJECTS_TSV, INVENTORY_TSV, AUDIT_TSV, LOG_DIR, WORK
EOF
  exit 1 ;;
esac
