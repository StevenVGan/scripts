#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0_config.sh — csRNA / nascent RNA fork (stranded tracks, optional post-trim MultiQC)
#
# Purpose
#   Central configuration file for the csRNA upstream pipeline (based on CUT&RUN cutrun).
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
#   - By default this file prepends conda env "bio" to PATH (see CONDA_BIO_ENV below)
#     so trim_galore, cutadapt, fastqc, bowtie2, samtools, deeptools (bamCoverage),
#     homer (makeTagDirectory), bedtools, macs3, etc. resolve without `conda activate`.
#   - To skip: CONDA_BIO_ENV="" ./run_all.sh  — or set CONDA_BIO_ENV to another env root.
#
# Notes
#   - This script creates the required output directories if they do not exist.
#   - Edit this file first when moving the project or changing genome/index paths.
###############################################################################


# ---- Base ----
# Default: combined HepG2 240216 + 240701 csRNA project (override with BASE=... ./run_all.sh)
BASE="${BASE:-${HOME}/work/seq/PROseq/240216_240701_csRNA_HepG2-ER3XHA_Steven}"

# ---- Conda env "bio" (prepend to PATH when directory exists) ----
CONDA_BIO_ENV="${CONDA_BIO_ENV:-${HOME}/miniforge3/envs/bio}"
if [[ -n "${CONDA_BIO_ENV}" && -d "${CONDA_BIO_ENV}/bin" ]]; then
  export PATH="${CONDA_BIO_ENV}/bin:${PATH}"
fi

# ---- Directories ----
RAW_DIR="${BASE}/data"              # 1_trim_qc input
TRIM_DIR="${BASE}/cleandata"        # 1_trim_qc out; 2_bowtie2 input
ALIGN_DIR="${BASE}/align"
BAM_DIR="${ALIGN_DIR}/bam"           # 2, 3, 4.1, 4.2, 5_qc
TRACK_DIR="${ALIGN_DIR}/track"       # 2_bowtie2 bigWigs
BAMQC_DIR="${BAM_DIR}/bamqc"        # 2 (idxstats, stats), 5_qc
TAG_DIR="${ALIGN_DIR}/tags"         # 3_homer_tags out; 4.2, 4.3 input
PEAK_DIR="${BASE}/peaks"
TSS_DIR="${BASE}/tss"                # 4.3_tss_csrna.sh out (findcsRNATSS.pl)
MACS3_DIR="${PEAK_DIR}/macs3"        # 4.1 out; 5_qc optional summary
HOMER_PEAK_DIR="${PEAK_DIR}/homer"   # 4.2 out
MULTIQC_OUT="${BASE}/multiqc"        # 5_qc
LOG_DIR="${BASE}/logs"

# ---- Single-end mode (SE=1) vs paired-end (SE=0, default) ----
# Affects: link_fastq (R2 optional), 1_trim_qc, 2_bowtie2, MACS3 format
SE="${SE:-0}"
[[ "$SE" -eq 1 ]] && MACS3_FORMAT="${MACS3_FORMAT:-BAM}" || MACS3_FORMAT="${MACS3_FORMAT:-BAMPE}"

# ---- Step toggles (used by run_all.sh) ----
# Preset via environment before ./run_all.sh overrides these defaults, e.g. RUN_TRIM=0 RUN_BOWTIE2=0 ./run_all.sh
RUN_TRIM="${RUN_TRIM:-1}"
RUN_BOWTIE2="${RUN_BOWTIE2:-1}"
RUN_HOMER_TAGS="${RUN_HOMER_TAGS:-1}"
RUN_PEAK_MACS3="${RUN_PEAK_MACS3:-0}"
RUN_PEAK_HOMER="${RUN_PEAK_HOMER:-1}"
RUN_TSS_CSRNA="${RUN_TSS_CSRNA:-1}" # 4.3_tss_csrna.sh: HOMER findcsRNATSS.pl — csRNA-specific TSS calling
RUN_QC="${RUN_QC:-1}"              # 5_qc.sh: preseq, plotFingerprint, PCA/Correlation, optional MACS3/TSS summary, MultiQC

# ---- Trimming (1_trim_qc) — csRNA / short 5′ RNAs ----
# Single Trim Galore call: auto-detects R1 3′ adapter, uses --adapter2 for R2's
# Illumina Small RNA 5′ readthrough, --nextseq for NovaSeq/NextSeq polyG tails.
# HOMER csRNA tutorial: reads shorter than ~15–20 nt are largely unusable; -min 20 fits 20–65 nt inserts.
TRIM_STRINGENCY=5
TRIM_MIN_LENGTH=20
# Quality cutoff passed via --nextseq (not -q). On NovaSeq/NextSeq 2-color chemistry, "no signal"
# is called as high-quality G; --nextseq ignores G qualities so polyG tails are trimmed like low quality.
TRIM_QUAL=25
TRIM_CPU=8
# Extra Trim Galore args (space-separated), e.g. "--clip_R1 2" — add after raw FastQC if needed
TRIM_GALORE_EXTRA="${TRIM_GALORE_EXTRA:-}"
# Illumina Small RNA 5′ adapter sequence — appears at R2 3′ end when insert < read length (csRNA
# inserts are 20–65 nt). Cutadapt post-step applies this as -A (R2-only). Set empty to skip.
SMALL_RNA_5P_ADAPTER="${SMALL_RNA_5P_ADAPTER:-GATCGTCGGACTGTAGAACTCTGAAC}"

# ---- Alignment (2_bowtie2) ----
BT2_CPU=12
BAMCOV_CPU=8
# Stranded PE orientation (required for correct strand bigWigs): set from kit chemistry, e.g. BT2_EXTRA="--fr" or BT2_EXTRA="--rf"
BT2_EXTRA="${BT2_EXTRA:-}"
# 1: emit strand-separated bigWigs (_fwd.bw / _rev.bw via bamCoverage --filterRNAstrand)
STRAND_BIGWIG="${STRAND_BIGWIG:-1}"
# When STRAND_BIGWIG=1: also emit unstranded ${sample}.bw (set 0 to save space)
COMBINED_BIGWIG="${COMBINED_BIGWIG:-1}"
# 1: HOMER makeTagDirectory -sspe (stranded PE); confirm against library prep
HOMER_SS_PE="${HOMER_SS_PE:-1}"
# ---- Reference (hg38) ----
GENOME="hg38"
GENOME_INDEX="/mnt/share/archive/bkup/ref/align/bowtie2/${GENOME}_noalt/${GENOME}"

# ---- QC (5_qc) ----
QC_THREADS=8
BAMQC_REGION="chr19"           # multiBamSummary: one chr for speed; empty = whole genome

# ---- Coverage / deepTools (2_bowtie2: bamCoverage) ----
BINSIZE=10
GENOMESIZE=2913022398          # deepTools hg38 effective genome size (matches GENOME=hg38)
NORMALIZE="RPGC"               # 1× depth scaling; requires --effectiveGenomeSize (deepTools docs)
IGNORE_CHR="chrM"
# HOMER csRNA guidance: do not routinely deduplicate csRNA (dynamic range like RNA-seq). 0 = do not pass --ignoreDuplicates to bamCoverage.
BAMCOV_IGNORE_DUP="${BAMCOV_IGNORE_DUP:-0}"

# ---- Space saving (2_bowtie2) ----
DELETE_TRIMMED_AFTER_ALIGN=1   # delete *_val_*.fq.gz after alignment

# ---- Peak calling (4.1 MACS3, 4.2 HOMER) ----
PEAKCALL_GROUPS_FILE="${BASE}/peakcall_groups.tsv"
# Format: ip_bam<TAB>control_bam<TAB>name<TAB>type  (control "-"/"none"/"NA" = no control; type TF|Histone)

BLACKLIST="${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed"

# MACS3 (4.1) — BAMPE for paired-end, BAM for single-end (set above from SE)
MACS3_CMD="macs3"
MACS3_GENOMESIZE="hs"          # MACS3 human preset (use with hg38/hg19 BAMs)
MACS3_FDR=0.05
MACS3_BROAD_CUTOFF=0.10

# ---- csRNA TSS calling (4.3_tss_csrna) ----
# 4.3 reads TSS_GROUPS_FILE — 3-col TSV (ip_bam<TAB>control_bam<TAB>name).
# Strict 1:1 mapping (no pooling / merging); rows with "-"/none/NA are skipped.
TSS_GROUPS_FILE="${TSS_GROUPS_FILE:-${BASE}/tss_groups.tsv}"
# HOMER findcsRNATSS.pl minimum tag threshold (Duttke et al. 2019 default = 7).
CSRNA_NTAG_THRESHOLD="${CSRNA_NTAG_THRESHOLD:-7}"
# Optional GTF for gene-annotation context (classification of promoter-proximal vs distal).
# Leave empty to rely on -genome $GENOME (HOMER built-in annotation).
CSRNA_GTF="${CSRNA_GTF:-}"

###############################################################################
# Helpers
###############################################################################
mkdir -p "$RAW_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$BAM_DIR" "$TRACK_DIR" "$BAMQC_DIR" "$TAG_DIR" "$PEAK_DIR" "$MACS3_DIR" "$HOMER_PEAK_DIR" "$TSS_DIR" "$LOG_DIR" "$MULTIQC_OUT"

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

  # Auto-prune: keep only LOG_KEEP_N most recent logs per step (default 3).
  # Logs under LOG_DIR/keep/ are exempt (glob doesn't recurse). See scripts/CONVENTIONS.md §10.
  local keep_n="${LOG_KEEP_N:-3}"
  if [[ -d "$LOG_DIR" && "$keep_n" =~ ^[0-9]+$ ]]; then
    ls -1t "${LOG_DIR}/${step}_"*.log 2>/dev/null | tail -n +$((keep_n + 1)) | xargs -r rm -f || true
  fi

  local logfile="${LOG_DIR}/${step}_$(date +%Y%m%d_%H%M%S).log"
  echo "Logging to: $logfile"
  exec > >(tee -a "$logfile") 2>&1
}

# Emit ${BASE}/references.tsv recording reference paths + env lock used.
# Called from 5_qc.sh as the last step. See scripts/CONVENTIONS.md §4.
emit_references_tsv() {
  local refs_file="${BASE}/references.tsv"
  local sz mt sha fasta lock_dir latest_lock

  {
    printf "key\tpath\tsize_bytes\tmtime_iso\tsha256_or_dash\n"

    printf "genome\t%s\t-\t-\t-\n" "${GENOME:-unknown}"

    fasta="${GENOME_INDEX}.fa"
    if [[ -f "$fasta" ]]; then
      sz=$(stat -c%s "$fasta" 2>/dev/null || echo "-")
      mt=$(date -u -d "@$(stat -c%Y "$fasta" 2>/dev/null)" -Iseconds 2>/dev/null || echo "-")
      printf "genome_fasta\t%s\t%s\t%s\t-\n" "$fasta" "$sz" "$mt"
    else
      printf "genome_fasta\t%s\t-\t-\t-\n" "$fasta"
    fi

    printf "bowtie2_index\t%s\t-\t-\t-\n" "$GENOME_INDEX"

    if [[ -f "$BLACKLIST" ]]; then
      sz=$(stat -c%s "$BLACKLIST" 2>/dev/null || echo "-")
      mt=$(date -u -d "@$(stat -c%Y "$BLACKLIST" 2>/dev/null)" -Iseconds 2>/dev/null || echo "-")
      sha=$(sha256sum "$BLACKLIST" 2>/dev/null | awk '{print $1}')
      [[ -z "$sha" ]] && sha="-"
      printf "blacklist\t%s\t%s\t%s\t%s\n" "$BLACKLIST" "$sz" "$mt" "$sha"
    else
      printf "blacklist\t%s\t-\t-\t-\n" "$BLACKLIST"
    fi

    lock_dir="${HOME}/work/scripts/env/lock"
    if [[ -d "$lock_dir" ]]; then
      latest_lock=$(ls -1t "$lock_dir"/bio.*.yml 2>/dev/null | head -1 || true)
      if [[ -n "$latest_lock" ]]; then
        sz=$(stat -c%s "$latest_lock" 2>/dev/null || echo "-")
        mt=$(date -u -d "@$(stat -c%Y "$latest_lock" 2>/dev/null)" -Iseconds 2>/dev/null || echo "-")
        printf "bio_env_lock\t%s\t%s\t%s\t-\n" "$latest_lock" "$sz" "$mt"
      fi
    fi
  } > "$refs_file"

  echo "references.tsv written: $refs_file"
}
