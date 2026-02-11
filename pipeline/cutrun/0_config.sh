#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0_config.sh
#
# Purpose
#   Central configuration file for the CUT&RUN upstream pipeline.
#   This file defines all project paths, parameters, and shared helper
#   functions used by the step scripts.
#
# What it defines
#   - Project directories (RAW/TRIM/ALIGN/BAM/TRACK/TAGS/PEAKS/LOGS)
#   - Step toggles used by run_all.sh (RUN_TRIM, RUN_BOWTIE2, RUN_HOMER_TAGS, RUN_PEAK_MACS3, RUN_PEAK_HOMER, RUN_QC)
#   - Trimming parameters (Trim Galore / cutadapt)
#   - Alignment parameters (bowtie2 index, CPU threads)
#   - Track/QC parameters (bamCoverage)
#   - MACS3 parameters (groups file, q-value, broad/narrow behavior, blacklist)
#   - Helper functions:
#       * check_cmd            : verify tool exists on PATH
#       * check_cmd_string     : verify first token for compound commands (e.g. "conda run ...")
#       * log_start            : create a timestamped log file and tee stdout/stderr into it
#
# Expected environment
#   - Run inside the conda env where core tools are available on PATH:
#       trim_galore, cutadapt, fastqc, bowtie2, samtools, deeptools (bamCoverage),
#       homer (makeTagDirectory), bedtools, and macs3 (or conda env "macs3").
#
# Notes
#   - This script creates the required output directories if they do not exist.
#   - Edit this file first when moving the project or changing genome/index paths.
###############################################################################


# ---- Base ----
# Set BASE for your project, or override via env: BASE=/path/to/project ./run_all.sh
BASE="${BASE:-${HOME}/work/seq/CUTRUN/my_project}"

# ---- Directories ----
RAW_DIR="${BASE}/data"              # 1_trim_qc input
TRIM_DIR="${BASE}/cleandata"        # 1_trim_qc out; 2_bowtie2 input
ALIGN_DIR="${BASE}/align"
BAM_DIR="${ALIGN_DIR}/bam"           # 2, 3, 4.1, 4.2, 5_qc
TRACK_DIR="${ALIGN_DIR}/track"       # 2_bowtie2 bigWigs
BAMQC_DIR="${BAM_DIR}/bamqc"        # 2 (idxstats, stats), 5_qc
TAG_DIR="${ALIGN_DIR}/tags"         # 3_homer_tags out; 4.2 input
PEAK_DIR="${BASE}/peaks"
MACS3_DIR="${PEAK_DIR}/macs3"        # 4.1 out; 5_qc optional summary
HOMER_PEAK_DIR="${PEAK_DIR}/homer"   # 4.2 out
MULTIQC_OUT="${BASE}/multiqc"        # 5_qc
LOG_DIR="${BASE}/logs"

# ---- Step toggles (used by run_all.sh) ----
RUN_TRIM=1
RUN_BOWTIE2=1
RUN_HOMER_TAGS=1
RUN_PEAK_MACS3=0
RUN_PEAK_HOMER=1
RUN_QC=1              # 5_qc.sh: preseq, run_spp.R, plotFingerprint, PCA/Correlation, optional MACS3 summary, MultiQC


# ---- Trimming (1_trim_qc) ----
TRIM_STRINGENCY=5
TRIM_MIN_LENGTH=36
TRIM_QUAL=25
TRIM_CPU=8

# ---- Alignment (2_bowtie2) ----
BT2_CPU=12
BAMCOV_CPU=8
GENOME="hg38"
GENOME_INDEX="/mnt/share/archive/bkup/ref/align/bowtie2/${GENOME}_noalt/${GENOME}"

# ---- QC (5_qc) ----
QC_THREADS=8
BAMQC_REGION="chr19"           # multiBamSummary: one chr for speed; empty = whole genome

# ---- Coverage / deepTools (2_bowtie2: bamCoverage) ----
BINSIZE=10
GENOMESIZE=2913022398          # hg38 effective genome size
NORMALIZE="RPGC"
IGNORE_CHR="chrM"

# ---- Space saving (2_bowtie2) ----
DELETE_TRIMMED_AFTER_ALIGN=1   # delete *_val_*.fq.gz after alignment

# ---- Peak calling (4.1 MACS3, 4.2 HOMER) ----
PEAKCALL_GROUPS_FILE="${BASE}/peakcall_groups.tsv"
# Format: ip_bam<TAB>control_bam<TAB>name<TAB>type  (control "-"/"none"/"NA" = no control; type TF|Histone)

BLACKLIST="${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed"

# MACS3 (4.1)
MACS3_CMD="macs3"
MACS3_GENOMESIZE="hs"
MACS3_FORMAT="BAMPE"
MACS3_FDR=0.05
MACS3_BROAD_CUTOFF=0.10

###############################################################################
# Helpers
###############################################################################
mkdir -p "$TRIM_DIR" "$ALIGN_DIR" "$BAM_DIR" "$TRACK_DIR" "$BAMQC_DIR" "$TAG_DIR" "$PEAK_DIR" "$MACS3_DIR" "$HOMER_PEAK_DIR" "$LOG_DIR" "$MULTIQC_OUT"

check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
}

# For commands like: "conda run -n macs3 macs3" (not a single binary)
check_cmd_string() {
  local cmdstr="$1"
  # Check first token only (e.g., "conda")
  local first
  first="$(echo "$cmdstr" | awk '{print $1}')"
  if ! command -v "$first" >/dev/null 2>&1; then
    echo "ERROR: Required command '$first' not found in PATH (needed for: $cmdstr)" >&2
    exit 1
  fi
}

log_start() {
  local step="$1"
  local logfile="${LOG_DIR}/${step}_$(date +%Y%m%d_%H%M%S).log"
  echo "Logging to: $logfile"
  exec > >(tee -a "$logfile") 2>&1
}
