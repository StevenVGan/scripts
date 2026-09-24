# Prep tools (upstream / project setup)

Scripts for downloading raw FASTQs, merging Illumina lanes, and linking into a CUTRUN (or other) project **`data/`** tree before **`1_trim_qc.sh`**.

| Script | Purpose |
|--------|---------|
| **download_igm_fastq_globus.sh** | **IGM runs (current):** pull a run IGM shared via Globus into `~/work/raw_seq/<RUN_ID>/` and `md5sum -c` it — `./download_igm_fastq_globus.sh <RUN_ID> "<link from the IGM share e-mail>"`; re-runnable (re-attaches / skips what is complete); `VERIFY_ONLY=1` after a pull made in the Globus web app; `BACKUP_ROOT=` copies the verified run to the archive volume |
| **start_gcp.sh** | start this node's Globus Connect Personal endpoint, detached and restricted to `~/work/raw_seq`; called by the script above (and by cron `@reboot`) |
| **download_fastq.sh** | LEGACY IGM FTP (`wget` + `md5sum -c`; IGM retired FTP in 2026); env-supply `FTP_USER`, `FTP_PASSWORD`, `FTP_BASE`, `DEST_DIR` |
| **link_fastq.sh** | Symlink from `RAW_DIR` using `prefix`→`newname` map (`*_R1_001` / `*_R2_001` or SE `.fastq.gz`) |
| **merge_lanes_inplace.sh** | Concatenate lanes in **`IGM_DIR`**, write `${sample}_R{1,2}.fastq.gz`, **delete** lane files; needs `IGM_DIR` + `MAP_FILE` |
| **link_merged_fastqs.sh** | Symlink merged `${sample}_R{1,2}.fastq.gz` into **`DEST_DIR`**; set **`MERGED_FASTQ_DIR`** or **`IGM_FASTQ_DIR`** (same folder); needs **`MAP_FILE`** |
| **download_geo_fastq_ena.sh** | ENA HTTPS `curl` for SRR list (SE-only per SRR); parallel **`DOWNLOAD_JOBS`** |

Typical **multi-lane IGM** flow:

```bash
# 1) Pull the run from IGM's Globus share (copy the link from IGM's e-mail; RUN_ID = the IGM run id):
./download_igm_fastq_globus.sh RUN_ID "https://app.globus.org/file-manager?origin_id=...&origin_path=..."

# 2) Merge + delete lanes (destructive in IGM_DIR)
IGM_DIR=~/work/raw_seq/RUN_ID MAP_FILE=~/work/seq/CUTRUN/proj/igm_to_sample.tsv ./merge_lanes_inplace.sh

# 3) Link into project data/
MERGED_FASTQ_DIR=~/work/raw_seq/RUN_ID DEST_DIR=~/work/seq/CUTRUN/proj/data \
  MAP_FILE=~/work/seq/CUTRUN/proj/igm_to_sample.tsv ./link_merged_fastqs.sh
```

**Single-lane** (or already one file per mate): skip merge; use **`link_fastq.sh`** with `link_sample.tsv`.

**GEO / ENA (SRR, single-end):** set `DEST_DIR` / `SRR_LIST_FILE` as env vars for **`download_geo_fastq_ena.sh`**, then run from this directory (or call with absolute path: `bash /path/to/prep/download_geo_fastq_ena.sh`).

## IGM via Globus

IGM delivers sequencing runs through Globus (FTP retired in 2026): the submitter's Globus identity (UCSD AD e-mail) or the lab's
share group gets a collection; you pull it onto this node. Nothing site-specific is stored in this repo — the collection id
comes from IGM's e-mail link each time, and the endpoint id is read from the local Globus Connect Personal install.

**One-time setup per node (≈ 20 min, two browser logins):**

```bash
# 1) globus-cli in its own env from the tracked yml, then log in once (prints a URL; paste the code back)
conda env create -f ~/work/scripts/env/globus.yml      # or: mamba create -n globus -c conda-forge globus-cli (latest, then re-export the yml)
~/miniforge3/envs/globus/bin/globus login && chmod 700 ~/.globus

# 2) Globus Connect Personal under $HOME (no root needed; the symlink is the stable GCP_HOME)
cd ~ && curl -O https://downloads.globus.org/globus-connect-personal/linux/stable/globusconnectpersonal-latest.tgz
dir=$(tar tzf globusconnectpersonal-latest.tgz | head -1 | cut -d/ -f1); tar xzf globusconnectpersonal-latest.tgz
ln -sfn "$dir" globusconnectpersonal && rm globusconnectpersonal-latest.tgz
mkdir -p -m 700 ~/.globusonline/lta && printf '~/work/raw_seq,0,1\n' > ~/.globusonline/lta/config-paths   # even a bare start stays restricted

# 3) register the endpoint (prints a URL; paste the code back), then lock the state dir down
~/globusconnectpersonal/globusconnectpersonal -setup -n "seq intake" --description "raw sequencing intake"
chmod -R go-rwx ~/.globusonline
~/miniforge3/envs/globus/bin/globus endpoint show "$(~/miniforge3/envs/globus/bin/globus endpoint local-id)" | grep Visibility   # must print "Visibility: False" (= not public)

# 4) start it (idempotent; also add to crontab:  @reboot sleep 120 && ~/work/scripts/pipeline/tools/prep/start_gcp.sh)
~/work/scripts/pipeline/tools/prep/start_gcp.sh
```

If the IGM collection is High Assurance (`globus collection show <uuid>` → "High Assurance: True"), register with
`--high-assurance --atm <minutes> --owner <your Globus identity>` instead, and expect `globus session update` when a session times out.

**Per run:** `./download_igm_fastq_globus.sh <RUN_ID> "<link>"` — or, when IGM shares one stable lab collection laid out as
`/<lab>/<RUN_ID>/`, put `export IGM_SOURCE_ROOT="<COLLECTION_UUID>:/<lab>"` in your `~/.bashrc` (never in this repo) and run just
`./download_igm_fastq_globus.sh <RUN_ID>`. It starts the endpoint if needed, lists the source, submits the
transfer with checksum verification, waits (Ctrl-C is safe — re-run later to re-attach), hard-fails on an `md5sum.txt` mismatch,
and writes `logs/download_complete.tsv` inside the run dir. A mapped collection asks once for a consent: run the printed
`globus session consent ...` command and re-run. Long-read / BAM deliveries: same command; expect the "no md5sum.txt" warning.
Web-app alternative: transfer the folder in app.globus.org onto the endpoint at `/~/work/raw_seq/<RUN_ID>/`, then
`VERIFY_ONLY=1 ./download_igm_fastq_globus.sh <RUN_ID>`.

**Known delivery quirks (2026-09):** the lab's share is one guest collection laid out as `/<lab>/<RUN_ID>/` with
`md5sum.txt`, the FASTQs, `Undetermined_*` and BCL Convert's `Logs/` + `Reports/`; `Reports/legacy/` is unreadable
(mode 0750 on IGM's side), so the script transfers with `--skip-source-errors` and reports skipped entries (a skipped
FASTQ is an error, anything else a warning). If Globus refuses a submit with "identical paths has not yet completed",
the script re-attaches to that task by its label. Throughput observed: ~100 MB/s. The manifest check deliberately avoids `md5sum --ignore-missing` (coreutils 8.25 silently drops an existing file with it).

**Endpoint status / stop:** `~/globusconnectpersonal/globusconnectpersonal -status` / `-stop`. **Upgrade:** download the latest
tarball, untar, `-stop`, re-point the `~/globusconnectpersonal` symlink, `start_gcp.sh` (`~/.globusonline` is untouched).

**Before pushing this repo:** `setup/audit_public_tree.sh` (refuses UUIDs, institutional/gmail e-mail addresses and secret literals in the lines a push would add, commit messages included).

**Legacy FTP (`download_fastq.sh`):** kept for reference; all credentials and the run id are env-supplied.

**Paths:** All scripts in this folder are path-agnostic; use absolute paths in **`MAP_FILE`** / env vars if you do not `cd` here.

Parent index: [../README.md](../README.md).
