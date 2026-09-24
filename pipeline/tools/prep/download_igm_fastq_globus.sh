#!/usr/bin/env bash
# Download one IGM sequencing run shared via Globus into $HOME/work/raw_seq/<RUN_ID>/ and verify it.
# Location: pipeline/tools/prep/download_igm_fastq_globus.sh
#
# IGM delivers by Globus (FTP retired in 2026). One-time setup per node (see README.md "IGM via Globus"):
# Globus Connect Personal registered (start_gcp.sh), globus-cli in the `globus` env, `globus login` done once.
#
# Usage:  ./download_igm_fastq_globus.sh <RUN_ID> "<SOURCE>"
#   RUN_ID  destination dir name under DEST_ROOT (the IGM run id, e.g. YYMMDD_<instrument>_<runnum>_<flowcell>)
#   SOURCE  where the run lives on the IGM collection, in any of these forms:
#             - the link from IGM's share e-mail (app.globus.org/file-manager?origin_id=<uuid>&origin_path=<path>)
#             - <COLLECTION_UUID>:/path/to/run/
#             - <COLLECTION_UUID>              (path defaults to /)
#   Both may also be given as env vars (RUN_ID=... SOURCE=...), so a wrapper can export them and exec this.
#   If IGM shares one stable lab collection, export IGM_SOURCE_ROOT="<COLLECTION_UUID>:/<lab dir>" (e.g. in
#   ~/.bashrc — never in this repo) and SOURCE may be omitted: it becomes $IGM_SOURCE_ROOT/<RUN_ID>/.
#
# Re-runnable: a re-run re-attaches to a transfer still in flight, reports a run that is already complete,
# and re-downloads only files that changed (checksum sync). Ctrl-C or a dropped SSH session never cancels
# the Globus task — re-run the same command later to finish the verification.
#
# Overrides: DEST_ROOT WAIT STALL_LOOPS NOTIFY VERIFY_ONLY BACKUP_ROOT CONDA_GLOBUS_ENV GCP_HOME GCP_PATHS
#   VERIFY_ONLY=1  no transfer; only the md5 check + completion stamp (e.g. after a pull in the Globus web app)
#   WAIT=0         submit the transfer and exit; re-run later to wait and verify
#   BACKUP_ROOT=   if set, copy the verified run to $BACKUP_ROOT/<RUN_ID>/ (the lab's archive volume)
set -euo pipefail

# ==== CONFIG: defaults; override with env vars (no need to edit) ==============================
DEST_ROOT="${DEST_ROOT:-${HOME}/work/raw_seq}"                          # runs land in DEST_ROOT/<RUN_ID>/
CONDA_GLOBUS_ENV="${CONDA_GLOBUS_ENV:-${HOME}/miniforge3/envs/globus}"  # "" = use the globus CLI already on PATH
GCP_HOME="${GCP_HOME:-${HOME}/globusconnectpersonal}"                   # Globus Connect Personal install (symlink)
GCP_PATHS="${GCP_PATHS:-rw~/work/raw_seq}"                              # what the endpoint may touch (start_gcp.sh)
WAIT="${WAIT:-1}"                    # 1 = wait for the transfer here; 0 = submit and exit
STALL_LOOPS="${STALL_LOOPS:-18}"     # 10-min polls with no progress before giving up (18 = 3 h)
NOTIFY="${NOTIFY:-on}"               # Globus e-mail on completion: on|off|succeeded|failed|inactive
VERIFY_ONLY="${VERIFY_ONLY:-0}"      # 1 = md5 check + stamp only
BACKUP_ROOT="${BACKUP_ROOT:-}"       # optional archive copy target
IGM_SOURCE_ROOT="${IGM_SOURCE_ROOT:-}"  # optional "<COLLECTION_UUID>:/<lab dir>" → SOURCE defaults to $IGM_SOURCE_ROOT/<RUN_ID>/
# ==============================================================================================

RUN_ID="${1:-${RUN_ID:-}}"
SOURCE="${2:-${SOURCE:-}}"
PREP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UUID_RE='[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}'
LABEL="raw_seq ${RUN_ID//[^A-Za-z0-9_,-]/_}"   # Globus labels allow only letters, digits, space, - _ ,

die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "[WARN] $*" >&2; }
step() { echo "[$1/9] $2"; }
now()  { date -Iseconds; }

# ---- 1. run id, destination, logs; VERIFY_ONLY short-circuit -------------------------------------------------
[[ -n "$RUN_ID" ]] || die "usage: $(basename "$0") <RUN_ID> \"<SOURCE>\"   (or export RUN_ID/SOURCE)"
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || die "RUN_ID must match ^[A-Za-z0-9][A-Za-z0-9._-]*$ : '$RUN_ID'"
DEST_DIR="${DEST_ROOT}/${RUN_ID}"
case "$(realpath -m "$DEST_DIR")" in
  "$(realpath -m "$DEST_ROOT")"/*) ;;
  *) die "destination escapes DEST_ROOT: $DEST_DIR" ;;
esac
LOG_DIR="${DEST_DIR}/logs"
mkdir -p "$DEST_DIR" "$LOG_DIR" && chmod 700 "$LOG_DIR"
STAMP="${LOG_DIR}/download_complete.tsv"
INCOMPLETE="${LOG_DIR}/download_incomplete.tsv"
TASK_FILE="${LOG_DIR}/globus_task.id"
LIST_FILE="${LOG_DIR}/filelist.tsv"
XFER_LOG="${LOG_DIR}/globus_transfer.log"
exec > >(tee -a "${LOG_DIR}/download.log") 2>&1
# On kill <pid> / Ctrl-C: stop the child `globus task wait` (it runs in the background and is awaited with `wait`,
# so the signal is acted on immediately rather than after the child's own 10-min timeout) and release the lock.
# The Globus task itself is a server-side object and keeps running; re-run this command to re-attach.
trap 'pkill -P $$ 2>/dev/null || true; echo "interrupted — the Globus task keeps running; re-run this command to re-attach" >&2; exit 143' INT TERM
echo "== $(now)  $(basename "$0") RUN_ID=$RUN_ID"

listing_hash_on_disk() {   # what the stamp recorded (column 6), "-" for web-app pulls, "" if no stamp
  [[ -f "$STAMP" ]] && awk -F'\t' 'NR==2{print $6}' "$STAMP" || true
}

archive_copy() {   # rsync the verified run to $BACKUP_ROOT/<RUN_ID>/ and md5-check the copy; non-zero on any failure
  if ! mkdir -p "${BACKUP_ROOT}/${RUN_ID}" 2>/dev/null || [[ ! -w "${BACKUP_ROOT}/${RUN_ID}" ]]; then
    echo "ERROR: BACKUP_ROOT is not writable: ${BACKUP_ROOT}/${RUN_ID}" >&2; return 1
  fi
  rsync -a --exclude 'logs/' "${DEST_DIR}/" "${BACKUP_ROOT}/${RUN_ID}/" || { echo "ERROR: rsync to ${BACKUP_ROOT}/${RUN_ID}/ failed" >&2; return 1; }
  if [[ -f "${DEST_DIR}/md5sum.txt" ]]; then
    (cd "${BACKUP_ROOT}/${RUN_ID}" && md5sum -c --quiet md5sum.txt) > "${LOG_DIR}/backup_md5check.log" 2>&1 \
      || { cat "${LOG_DIR}/backup_md5check.log" >&2; echo "ERROR: archive copy failed its md5 check" >&2; return 1; }
  fi
  echo "archive copy verified: ${BACKUP_ROOT}/${RUN_ID}/"
}

verify_and_stamp() {       # step 8 (+9): md5 hard-fail, completion stamp, optional archive copy
  local task="${1:-$( [[ -f "$TASK_FILE" ]] && cat "$TASK_FILE" || echo '-')}"
  local list_hash="${2:--}" md5="no-manifest" n bytes backup="-"
  step 8 "verify"
  if [[ -f "${DEST_DIR}/md5sum.txt" ]]; then
    # Every manifest entry must be on disk (checked explicitly — coreutils 8.25's `md5sum --ignore-missing`
    # can silently drop an existing file), then a plain md5sum -c over the manifest.
    local listed missing
    listed=$(grep -cvE '^\s*(#|$)' "${DEST_DIR}/md5sum.txt" || true)
    missing=$(grep -vE '^\s*(#|$)' "${DEST_DIR}/md5sum.txt" | sed -E 's/^[0-9A-Fa-f]{32}[ *]+//' \
              | while IFS= read -r f; do [[ -e "${DEST_DIR}/${f}" ]] || printf '%s\n' "$f"; done)
    if [[ -n "$missing" ]]; then
      printf '%s\n' "$missing" > "${LOG_DIR}/md5check.log"
      echo "ERROR: $(grep -c '' <<<"$missing") of the $listed files in md5sum.txt are not on disk (names in logs/md5check.log)" >&2
      return 1
    fi
    if (cd "$DEST_DIR" && md5sum -c md5sum.txt) > "${LOG_DIR}/md5check.log" 2>&1; then
      local ok; ok=$(grep -c ': OK$' "${LOG_DIR}/md5check.log" || true)
      md5="ok:${ok}/${listed}"
      (( ok == listed )) || { echo "ERROR: md5sum -c verified $ok of $listed manifest entries (see logs/md5check.log)" >&2; return 1; }
      echo "md5: $ok/$listed OK"
    else
      grep -v ': OK$' "${LOG_DIR}/md5check.log" | head -20 >&2
      return 1
    fi
  else
    warn "no md5sum.txt was delivered — integrity rests on Globus per-file checksum verification (on)"
  fi
  n=$(find "$DEST_DIR" -path "$LOG_DIR" -prune -o -type f -print | wc -l)
  bytes=$(find "$DEST_DIR" -path "$LOG_DIR" -prune -o -type f -printf '%s\n' | awk '{s+=$1} END{print s+0}')
  (( n > 0 )) || { echo "ERROR: nothing to stamp — $DEST_DIR holds no files" >&2; return 1; }
  if [[ ! -f "${DEST_DIR}/md5sum.txt" && -f "$LIST_FILE" ]]; then   # no manifest: every listed top-level file must be on disk
    local absent
    absent=$(awk -F'\t' '$1=="file"{print $2}' "$LIST_FILE" | while IFS= read -r f; do [[ -e "${DEST_DIR}/${f}" ]] || printf '%s\n' "$f"; done)
    [[ -z "$absent" ]] || { printf '%s\n' "$absent" > "${LOG_DIR}/md5check.log"; echo "ERROR: $(grep -c '' <<<"$absent") listed files are not on disk (names in logs/md5check.log)" >&2; return 1; }
  fi
  # The run is verified at this point: stamp it even if the archive copy below fails (re-run with VERIFY_ONLY=1
  # after fixing BACKUP_ROOT to retry the copy).
  local backup_rc=0
  if [[ -n "$BACKUP_ROOT" ]]; then
    step 9 "archive copy → ${BACKUP_ROOT}/${RUN_ID}/"
    if archive_copy; then backup="${BACKUP_ROOT}/${RUN_ID}"; else backup="FAILED"; backup_rc=1; fi
  else
    step 9 "no BACKUP_ROOT set"
    warn "single copy — this volume has no backup; set BACKUP_ROOT to copy the run to the archive volume"
  fi
  printf 'run_id\ttask_id\tcompleted_at\tn_files\tbytes\tlisting_sha256\tmd5\tbackup\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$RUN_ID" "$task" "$(now)" "$n" "$bytes" "$list_hash" "$md5" "$backup" > "${STAMP}.tmp" && mv -f "${STAMP}.tmp" "$STAMP"
  rm -f "$INCOMPLETE"
  echo "Done. $RUN_ID: $n files, $bytes bytes, md5 $md5, backup $backup, stamp $STAMP" >&2
  (( backup_rc == 0 )) || die "the run is verified and stamped, but the archive copy FAILED — fix BACKUP_ROOT, then re-run with VERIFY_ONLY=1 to retry the copy"
}


exec 9>"${LOG_DIR}/download.lock"
flock -n 9 || die "another run of this command holds ${LOG_DIR}/download.lock — tail ${LOG_DIR}/download.log, or wait for it"
if [[ "$VERIFY_ONLY" == "1" ]]; then
  step 1 "VERIFY_ONLY=1 → md5 check + stamp for $DEST_DIR"
  # keep the task id and listing hash of an earlier full run, if any, so a re-delivery can still be detected later
  prev_task="$( [[ -f "$STAMP" ]] && awk -F'\t' 'NR==2{print $2}' "$STAMP" || echo '-')"
  prev_hash="$(listing_hash_on_disk)"; [[ -n "$prev_hash" ]] || prev_hash="-"
  verify_and_stamp "${prev_task:--}" "$prev_hash" || die "verification failed for $DEST_DIR (see logs/md5check.log)"
  exit 0
fi
step 1 "destination $DEST_DIR"

# ---- 2. CLI, login, source, lock ---------------------------------------------------------------------------------
[[ -n "$CONDA_GLOBUS_ENV" && -d "${CONDA_GLOBUS_ENV}/bin" ]] && export PATH="${CONDA_GLOBUS_ENV}/bin:${PATH}"
command -v globus >/dev/null || die "globus CLI not found — create the env: mamba create -n globus -c conda-forge globus-cli"
set +e; globus whoami >/dev/null 2>&1; rc=$?; set -e
(( rc == 0 )) || { (( rc == 4 )) && die "not logged in to Globus — run:  globus login   (paste the code it asks for)"; die "globus whoami failed (rc=$rc): network or CLI problem"; }
[[ -n "$SOURCE" || -z "$IGM_SOURCE_ROOT" ]] || SOURCE="${IGM_SOURCE_ROOT%/}/${RUN_ID}/"
[[ -n "$SOURCE" ]] || die "SOURCE missing: the link from IGM's share e-mail, <COLLECTION_UUID>:/path/, or export IGM_SOURCE_ROOT"
if [[ "$SOURCE" == *origin_id=* ]]; then
  t="${SOURCE#*origin_id=}"; [[ "$t" =~ ^$UUID_RE ]] || die "origin_id in the link is not a UUID"
  SRC_ID="${BASH_REMATCH[0]}"
else
  [[ "$SOURCE" =~ $UUID_RE ]] || die "SOURCE contains no collection UUID: '$SOURCE'"
  SRC_ID="${BASH_REMATCH[0]}"
fi
if [[ "$SOURCE" == *origin_path=* ]]; then
  p="${SOURCE#*origin_path=}"; p="${p%%&*}"; SRC_PATH="$(printf '%b' "${p//%/\\x}")"
elif [[ "$SOURCE" == *"${SRC_ID}:"* ]]; then
  SRC_PATH="${SOURCE#*"${SRC_ID}:"}"
else
  SRC_PATH="/"
fi
[[ "$SRC_PATH" == /* ]] || SRC_PATH="/${SRC_PATH}"
[[ "$SRC_PATH" == */ ]] || SRC_PATH="${SRC_PATH}/"
SRC="${SRC_ID}:${SRC_PATH}"
step 2 "source $SRC"

# ---- 3. endpoint -------------------------------------------------------------------------------------------------
step 3 "Globus Connect Personal"
GCP_HOME="$GCP_HOME" GCP_PATHS="$GCP_PATHS" "${PREP_DIR}/start_gcp.sh"
LOCAL_ID="$(globus endpoint local-id)" || die "no local Globus Connect Personal endpoint id (~/.globusonline/lta/client-id.txt)"
DEST="${LOCAL_ID}:${DEST_DIR}/"
# the endpoint only accepts writes under GCP_PATHS (rw prefixes); refuse now rather than after 30 min of PERMISSION_DENIED
dest_allowed=0
IFS=',' read -ra gcp_rules <<<"$GCP_PATHS"
for rule in "${gcp_rules[@]}"; do
  [[ "$rule" =~ ^(rw|w)(.*)$ ]] || continue
  allowed="${BASH_REMATCH[2]}"; allowed="${allowed/#\~/$HOME}"
  case "$(realpath -m "$DEST_DIR")/" in "$(realpath -m "$allowed")"/*) dest_allowed=1 ;; esac
done
(( dest_allowed )) || die "DEST_DIR $DEST_DIR is outside the endpoint's writable paths (GCP_PATHS=$GCP_PATHS) — set GCP_PATHS to match DEST_ROOT and restart the endpoint"

# ---- 4. listing --------------------------------------------------------------------------------------------------
step 4 "listing $SRC"
consent_or_die() {   # $1 = rc, $2 = stderr file
  if (( $1 == 4 )); then
    cat "$2" >&2
    die "Globus needs a one-time consent for this collection — run the 'globus session consent ...' command printed above, then re-run"
  fi
  cat "$2" >&2; die "globus command failed (rc=$1)"
}
errf="$(mktemp)"; trap 'rm -f "$errf"' EXIT
set +e; globus ls "$SRC" --format unix --jmespath 'DATA[].[type, name, size]' > "${LIST_FILE}.new" 2>"$errf"; rc=$?; set -e
(( rc == 0 )) || consent_or_die "$rc" "$errf"
sort -o "${LIST_FILE}.new" "${LIST_FILE}.new"
n_files=$(awk -F'\t' '$1=="file"' "${LIST_FILE}.new" | wc -l)
n_dirs=$(awk -F'\t' '$1=="dir"' "${LIST_FILE}.new" | wc -l)
(( n_files > 0 )) || die "no files at $SRC ($n_dirs directories) — point origin_path at the run directory itself"
LIST_HASH="$(sha256sum "${LIST_FILE}.new" | cut -c1-64)"
mv -f "${LIST_FILE}.new" "$LIST_FILE"
echo "listing: $n_files files, $n_dirs directories, hash ${LIST_HASH:0:12}"

# ---- 5. what to do: already complete / re-attach / verify / submit ------------------------------------------------
submit() {
  step 6 "submit transfer → $DEST"
  set +e
  # --skip-source-errors: IGM run dirs carry unreadable sub-dirs (e.g. Reports/legacy/, mode 0750); without it the
  # task retries them for 3 days. Skipped entries are checked after completion (a skipped FASTQ is an error).
  TASK="$(globus transfer "$SRC" "$DEST" --recursive --sync-level checksum --verify-checksum --preserve-timestamp \
          --skip-source-errors --label "$LABEL" --notify "$NOTIFY" --format unix --jmespath task_id 2>"$errf")"; rc=$?
  set -e
  if (( rc != 0 )) && grep -q "identical paths has not yet completed" "$errf"; then
    # Globus refuses a second task with the same source/destination while one is in flight (e.g. the task id file was
    # lost, or the transfer was started from the web app). Find it by label and re-attach instead of failing.
    TASK="$(globus task list --filter-status ACTIVE --filter-label "$LABEL" --format unix --jmespath 'DATA[0].task_id' 2>/dev/null || true)"
    [[ "$TASK" =~ ^$UUID_RE$ ]] || { cat "$errf" >&2; die "a transfer with the same paths is still in flight but could not be found by label '$LABEL' — check app.globus.org/activity, then re-run"; }
    echo "an identical transfer is already in flight — re-attaching to task $TASK"
    echo "$TASK" > "$TASK_FILE"; return 0
  fi
  (( rc == 0 )) && [[ "$TASK" =~ ^$UUID_RE$ ]] || consent_or_die "$rc" "$errf"
  echo "$TASK" > "$TASK_FILE"; rm -f "$INCOMPLETE"
  echo "$(now)  submitted task $TASK  $SRC -> $DEST" >> "$XFER_LOG"
  echo "task $TASK"
}

task_field() { globus task show "$1" --format unix --jmespath "$2" 2>/dev/null || echo "?"; }

wait_task() {
  step 7 "waiting for task $TASK (10-min polls; Ctrl-C is safe — the task keeps running)"
  local prev="" stalls=0 faulting=0 transient=0 status nice bytes faults ftx fskip subok sig
  while :; do
    set +e
    globus task wait "$TASK" --timeout 600 --timeout-exit-code 50 --polling-interval 30 -H &   # background + wait: a signal interrupts at once
    wait $!; rc=$?
    set -e
    if (( rc == 0 )); then break; fi
    if (( rc != 50 )); then   # not a timeout: a transient CLI/API/network error, not a task state — retry, bounded
      (( ++transient <= 10 )) || die "globus task wait kept failing (rc=$rc, 10 times) — task $TASK left running; re-run to re-attach"
      warn "globus task wait failed (rc=$rc); retrying in 60 s ($transient/10)"; sleep 60; continue
    fi
    transient=0
    IFS=$'\t' read -r status nice bytes faults ftx fskip subok <<<"$(task_field "$TASK" '[status, nice_status, bytes_transferred, faults, files_transferred, files_skipped, subtasks_succeeded]')"
    echo "$(now)  $status nice_status=$nice bytes=$bytes files=$ftx skipped=$fskip subtasks_ok=$subok faults=$faults" >> "$XFER_LOG"
    sig="$bytes/$ftx/$fskip/$subok"   # a checksum-only pass moves no bytes but does advance these counters
    if [[ "$status" == "INACTIVE" ]]; then
      die "task $TASK is INACTIVE ($nice): credentials/consent lapsed — run 'globus login', the printed 'globus session consent ...', or 'globus session update' (High Assurance), then re-run"
    fi
    [[ "$nice" == "OK" || "$nice" == "None" || "$nice" == "null" || -z "$nice" ]] || warn "nice_status=$nice"
    if [[ "$nice" == PERMISSION_DENIED* || "$nice" == FILE_NOT_FOUND* ]]; then
      (( ++faulting > 2 )) && {
        globus task cancel "$TASK" >/dev/null 2>&1 || true
        printf 'run_id\ttask_id\tcancelled_at\treason\n%s\t%s\t%s\t%s\n' "$RUN_ID" "$TASK" "$(now)" "$nice" > "$INCOMPLETE"
        die "task $TASK kept failing with $nice for >30 min — cancelled. The share may be gone (contact IGM) or the endpoint's path restriction blocks the write"
      }
    else
      faulting=0
    fi
    if [[ "$sig" == "$prev" ]]; then
      (( ++stalls ))
      (( stalls >= 6 )) && { warn "no progress for $stalls polls"; GCP_HOME="$GCP_HOME" GCP_PATHS="$GCP_PATHS" "${PREP_DIR}/start_gcp.sh" || true; }
      (( stalls >= STALL_LOOPS )) && die "transfer stalled for $stalls polls; task $TASK left running — re-run to re-attach"
    else
      stalls=0; prev="$sig"
    fi
  done
  status="$(task_field "$TASK" status)"
  { echo "$(now)  final $status"; globus task show "$TASK" 2>/dev/null || true; } >> "$XFER_LOG"
  if [[ "$status" != "SUCCEEDED" ]]; then
    globus task event-list "$TASK" --filter-errors --limit 20 >> "$XFER_LOG" 2>&1 || true
    die "task $TASK ended with status $status (details: $XFER_LOG)"
  fi
  check_skipped "$TASK"
  echo "task $TASK SUCCEEDED"
}

check_skipped() {   # entries Globus skipped on the source (--skip-source-errors): a FASTQ is fatal, anything else a warning
  local skipped
  skipped="$(globus task show "$1" --skipped-errors --format unix --jmespath 'DATA[].[error_code, source_path]' 2>/dev/null || true)"
  [[ -n "$skipped" && "$skipped" != "None" ]] || return 0
  printf 'skipped source entries (task %s):\n%s\n' "$1" "$skipped" >> "$XFER_LOG"
  local files_skipped
  files_skipped="$(awk -F'\t' '$2 !~ /\/$/ {print $2}' <<<"$skipped")"   # directories end with "/", files do not
  if [[ -n "$files_skipped" ]]; then
    printf '%s\n' "$files_skipped" >&2; die "a data FILE on the source was skipped (unreadable) — contact IGM"
  fi
  warn "skipped unreadable source directories (see $XFER_LOG): $(awk -F'\t' '{print $2}' <<<"$skipped" | tr '\n' ' ')"
}

verify_or_repair() {   # after a SUCCEEDED task: verify; on md5 mismatch re-sync once (checksum sync repairs bad files)
  if verify_and_stamp "$TASK" "$LIST_HASH"; then return 0; fi
  (( ${REPAIRED:-0} == 0 )) || die "md5 still failing after a re-transfer — the manifest or the source is wrong; contact IGM"
  warn "md5 mismatch — re-syncing with checksum comparison to repair the mismatching file(s)"
  REPAIRED=1; submit; wait_task; verify_or_repair
}

step 5 "state"
prev_hash="$(listing_hash_on_disk)"
if [[ -n "$prev_hash" ]]; then
  if [[ "$prev_hash" == "$LIST_HASH" ]]; then
    echo "already complete: stamp $STAMP matches the current listing"; exit 0
  elif [[ "$prev_hash" == "-" ]]; then
    missing=$(awk -F'\t' '$1=="file"{print $2}' "$LIST_FILE" | while read -r f; do [[ -e "${DEST_DIR}/${f}" ]] || echo "$f"; done | wc -l)
    (( missing == 0 )) && { echo "already complete: every listed file is on disk (web-app pull, stamped)"; exit 0; }
    echo "$missing listed files are missing on disk — transferring"
  else
    echo "the listing changed since the stamp (re-delivery?) — transferring what is new"
  fi
fi
if [[ -f "$TASK_FILE" ]]; then
  TASK="$(cat "$TASK_FILE")"
  status="$(task_field "$TASK" status)"
  case "$status" in
    ACTIVE)    echo "re-attaching to task $TASK (ACTIVE)"; [[ "$WAIT" == "1" ]] || { echo "WAIT=0: task still running; re-run later"; exit 0; }
               wait_task; verify_or_repair; exit 0 ;;
    INACTIVE)  die "task $TASK is INACTIVE ($(task_field "$TASK" nice_status)): run 'globus login' / the printed 'globus session consent ...' / 'globus session update', then re-run" ;;
    SUCCEEDED) [[ -n "$prev_hash" && "$prev_hash" != "-" ]] || { echo "task $TASK SUCCEEDED earlier — verifying"; check_skipped "$TASK"; verify_or_repair; exit 0; } ;;
    *)         echo "recorded task $TASK is $status — submitting a new transfer" ;;
  esac
fi
submit
[[ "$WAIT" == "1" ]] || { echo "WAIT=0: submitted; re-run the same command later to wait and verify"; exit 0; }
wait_task
verify_or_repair
