#!/usr/bin/env bash
# Start this node's Globus Connect Personal endpoint, detached and restricted to GCP_PATHS. Idempotent.
# Location: pipeline/tools/prep/start_gcp.sh
# Usage:  ./start_gcp.sh          (download_igm_fastq_globus.sh calls it; also from cron:  @reboot sleep 120 && ~/work/scripts/pipeline/tools/prep/start_gcp.sh)
#   status:  ~/globusconnectpersonal/globusconnectpersonal -status
#   stop:    ~/globusconnectpersonal/globusconnectpersonal -stop
# Overrides: GCP_HOME GCP_PATHS GCP_WAIT_SECS
set -euo pipefail

# ==== CONFIG ==================================================================================
GCP_HOME="${GCP_HOME:-${HOME}/globusconnectpersonal}"   # symlink to the unpacked globusconnectpersonal-<version>/
GCP_PATHS="${GCP_PATHS:-rw~/work/raw_seq}"             # -restrict-paths: the ONLY paths Globus may touch on this node
GCP_WAIT_SECS="${GCP_WAIT_SECS:-90}"                   # how long to wait for "connected"
# ==============================================================================================

BIN="${GCP_HOME}/globusconnectpersonal"
STATE="${HOME}/.globusonline"
die() { echo "ERROR: $*" >&2; exit 1; }

[[ -x "$BIN" ]] || die "Globus Connect Personal not found at $BIN (see prep/README.md 'IGM via Globus')"
for f in config client-id.txt gridmap; do
  [[ -f "${STATE}/lta/${f}" ]] || die "endpoint not registered yet — run:  $BIN -setup -n \"<endpoint name>\" --description \"<text>\""
done

gcp_connected() { local out; out="$(timeout 20 "$BIN" -status 2>/dev/null || true)"; grep -Eq '^Globus Online:[[:space:]]+connected' <<<"$out"; }
gcp_running()   { timeout 20 "$BIN" -status >/dev/null 2>&1; }   # exit 0 whenever an instance holds the control socket

if gcp_connected; then echo "[gcp] running and connected"; exit 0; fi
if gcp_running; then
  echo "[gcp] running but not connected — waiting up to ${GCP_WAIT_SECS}s"
  for ((i = 0; i < GCP_WAIT_SECS; i += 5)); do sleep 5; gcp_connected && { echo "[gcp] connected"; exit 0; }; done
  echo "[gcp] still not connected — restarting"
  "$BIN" -stop >/dev/null 2>&1 || true; sleep 3
fi

# Clean environment: GCP's own launcher picks `python3` from PATH (system python is fine, a conda env is not wanted).
# Close every inherited descriptor above 2 first: when this is called from inside download_igm_fastq_globus.sh the
# daemon would otherwise keep that run's lock file and log open for its whole lifetime.
(
  for fdpath in /proc/$BASHPID/fd/*; do n="${fdpath##*/}"; (( n > 2 && n != 255 )) && eval "exec ${n}>&-"; done
  echo "== $(date -Iseconds) start_gcp.sh: starting $BIN -start -restrict-paths $GCP_PATHS" >> "${STATE}/gcp.nohup"
  exec env -i HOME="$HOME" USER="${USER:-$(id -un)}" LOGNAME="${LOGNAME:-$(id -un)}" PATH=/usr/bin:/bin \
    setsid nohup "$BIN" -start -restrict-paths "$GCP_PATHS" >> "${STATE}/gcp.nohup" 2>&1 < /dev/null
) &
for ((i = 0; i < GCP_WAIT_SECS; i += 5)); do
  sleep 5
  gcp_connected && { echo "[gcp] started and connected (restricted to: $GCP_PATHS)"; exit 0; }
done
tail -20 "${STATE}/gcp.nohup" >&2 || true
die "Globus Connect Personal did not reach 'connected' within ${GCP_WAIT_SECS}s (see ${STATE}/gcp.nohup)"
