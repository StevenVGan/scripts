#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0_config.sh — PRO-seq / GRO-seq fork (stranded tracks + strand flip + pausing index)
#
# Purpose
#   Central configuration file for the PRO-seq upstream pipeline.
#   Forked from scripts/pipeline/csRNA/ with three PRO-seq-specific changes:
#     1. PROSEQ_FLIP_STRAND=1 — PRO-seq R1 is the reverse complement of the
#        nascent RNA (Kwak et al. 2013), so tag directories and stranded bigWigs
#        must flip strand to represent the RNA strand.
#     2. HOMER_STYLE=groseq — findPeaks -style groseq replaces the csRNA-only
#        findcsRNATSS.pl; appropriate for nascent transcript detection.
#     3. Step 4.3 replaced by pausing index + divergent eRNA analysis
#        (scripts/pipeline/proseq/4.3_pausing_divergent.sh).
#
# Expected environment
#   - Prepends conda env "bio" to PATH when available (trim_galore, cutadapt,
#     fastqc, bowtie2, samtools, deeptools, homer, bedtools, macs3).
#   - Skip: CONDA_BIO_ENV="" ./run_all.sh
###############################################################################


# ---- Base ----
BASE="${BASE:-${HOME}/work/seq/PROseq/my_proseq_project}"

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
BAMQC_DIR="${BAM_DIR}/bamqc"        # 2, 5_qc
TAG_DIR="${ALIGN_DIR}/tags"         # 3_homer_tags out; 4.2, 4.3 input
PEAK_DIR="${BASE}/peaks"
MACS3_DIR="${PEAK_DIR}/macs3"        # 4.1 out
HOMER_PEAK_DIR="${PEAK_DIR}/homer"   # 4.2 out (findPeaks -style groseq)
PAUSING_DIR="${BASE}/pausing"        # 4.3 out (pausing index + divergent pairs)
MULTIQC_OUT="${BASE}/multiqc"        # 5_qc
LOG_DIR="${BASE}/logs"

# ---- Single-end mode (SE=1) vs paired-end (SE=0, default) ----
SE="${SE:-0}"
[[ "$SE" -eq 1 ]] && MACS3_FORMAT="${MACS3_FORMAT:-BAM}" || MACS3_FORMAT="${MACS3_FORMAT:-BAMPE}"

# ---- Step toggles (used by run_all.sh) ----
RUN_TRIM="${RUN_TRIM:-1}"
RUN_POLYA_TRIM="${RUN_POLYA_TRIM:-0}"  # 1.1_polya_trim.sh — enable per project when FastQC flags poly-A/T
RUN_BOWTIE2="${RUN_BOWTIE2:-1}"
RUN_HOMER_TAGS="${RUN_HOMER_TAGS:-1}"
RUN_PEAK_MACS3="${RUN_PEAK_MACS3:-0}"  # rarely appropriate for nascent RNA
RUN_PEAK_HOMER="${RUN_PEAK_HOMER:-1}"  # findPeaks -style groseq (transcript calls)
RUN_PAUSING="${RUN_PAUSING:-1}"        # 4.3_pausing_divergent.sh
RUN_QC="${RUN_QC:-1}"

# ---- Poly-A / poly-T trim (1.1) ------------------------------------------
# Only consulted when RUN_POLYA_TRIM=1. Defaults match cutadapt's usual.
POLYA_MIN_STRETCH="${POLYA_MIN_STRETCH:-10}"    # minimum A/T run length to trim
POLYA_ERROR_RATE="${POLYA_ERROR_RATE:-0.1}"     # cutadapt -e

# ---- PRO-seq strand handling ---------------------------------------------
# PRO-seq R1 is the reverse complement of the nascent RNA (biotin-NTP run-on →
# streptavidin pulldown → small-RNA library prep places the 5′ Illumina adapter
# on the RNA's 3′ end). R1 therefore maps to the OPPOSITE strand of the RNA.
#   1 → HOMER makeTagDirectory gets -flip, and bamCoverage --filterRNAstrand
#       labels are SWAPPED so ${sample}_fwd.bw = RNA on + strand.
#   0 → leave strand as-aligned (csRNA-style; useful for PRO-cap variants).
PROSEQ_FLIP_STRAND="${PROSEQ_FLIP_STRAND:-1}"

# ---- Trimming (1_trim_qc) — PRO-seq small-RNA chemistry ------------------
# PRO-seq inserts ~30–40 nt (similar to csRNA); min length 20 avoids garbage.
# --nextseq handles polyG tails from NovaSeq/NextSeq 2-color chemistry.
TRIM_STRINGENCY=5
TRIM_MIN_LENGTH=20
TRIM_QUAL=25
TRIM_CPU=8
TRIM_GALORE_EXTRA="${TRIM_GALORE_EXTRA:-}"
# Illumina Small RNA 5′ adapter — appears at R2 3′ end when insert < read length.
# PRO-seq uses small-RNA chemistry so the same readthrough concern applies.
SMALL_RNA_5P_ADAPTER="${SMALL_RNA_5P_ADAPTER:-GATCGTCGGACTGTAGAACTCTGAAC}"

# ---- Alignment (2_bowtie2) -----------------------------------------------
BT2_CPU=12
BAMCOV_CPU=8
# Stranded PE orientation from library chemistry: TruSeq-style PRO-seq = "--fr".
BT2_EXTRA="${BT2_EXTRA:-}"
# 1: emit strand-separated bigWigs (_fwd.bw / _rev.bw via bamCoverage --filterRNAstrand).
STRAND_BIGWIG="${STRAND_BIGWIG:-1}"
# When STRAND_BIGWIG=1: also emit unstranded ${sample}.bw.
COMBINED_BIGWIG="${COMBINED_BIGWIG:-1}"
# 1: HOMER makeTagDirectory -sspe (stranded PE). Forced off for SE.
HOMER_SS_PE="${HOMER_SS_PE:-1}"
[[ "$SE" -eq 1 ]] && HOMER_SS_PE=0

# ---- Reference (hg38) ----
GENOME="${GENOME:-hg38}"
GENOME_INDEX="/mnt/share/archive/bkup/ref/align/bowtie2/${GENOME}_noalt/${GENOME}"

# ---- QC (5_qc) ----
QC_THREADS=8
BAMQC_REGION="chr19"

# ---- Coverage / deepTools (2_bowtie2: bamCoverage) -----------------------
BINSIZE=10
GENOMESIZE=2913022398          # deepTools hg38 effective genome size
NORMALIZE="RPGC"
IGNORE_CHR="chrM"
# PRO-seq dynamic range resembles RNA-seq — do not routinely deduplicate.
BAMCOV_IGNORE_DUP="${BAMCOV_IGNORE_DUP:-0}"

# ---- Space saving (2_bowtie2) --------------------------------------------
DELETE_TRIMMED_AFTER_ALIGN=1

# ---- Peak calling (4.1 MACS3, 4.2 HOMER) ---------------------------------
# Groups file format: ip_bam<TAB>control_bam<TAB>name<TAB>type
#   (control "-"/"none"/"NA" = no control; type TF|Histone — ignored by -style groseq)
PEAKCALL_GROUPS_FILE="${PEAKCALL_GROUPS_FILE:-${BASE}/peakcall_groups.tsv}"
BLACKLIST="${BLACKLIST:-${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed}"

# MACS3 (4.1) — off by default; kept available for specific analyses.
MACS3_CMD="macs3"
MACS3_GENOMESIZE="hs"
MACS3_FDR=0.05
MACS3_BROAD_CUTOFF=0.10

# HOMER findPeaks (4.2) — groseq style for PRO-seq nascent transcript calling.
# Other valid styles: factor, histone, tss, mRNA, super, etc.
HOMER_STYLE="${HOMER_STYLE:-groseq}"

# ---- Pausing index + divergent TSS (4.3) ---------------------------------
# Annotation: BED12 transcript file (Steven's convention at ~/work/ref/annot/<genome>/).
GENE_ANNOT_BED="${GENE_ANNOT_BED:-${HOME}/work/ref/annot/${GENOME}/${GENOME}_annot_genome.bed}"
# Promoter window half-width (bp around TSS) for the numerator of the pausing index.
TSS_WINDOW="${TSS_WINDOW:-250}"
# Gene body starts this many bp downstream of TSS (excluding promoter region).
GB_START="${GB_START:-500}"
# Minimum gene length (bp) to include (need meaningful gene body downstream of GB_START).
GB_MIN_LEN="${GB_MIN_LEN:-1000}"
# Divergent eRNA detection: max distance (bp) upstream of a + transcript's 5′
# end to accept an antisense (- strand) call from findPeaks as its divergent pair.
DIVERGENT_WINDOW="${DIVERGENT_WINDOW:-1000}"

###############################################################################
# Helpers
###############################################################################
mkdir -p "$RAW_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$BAM_DIR" "$TRACK_DIR" "$BAMQC_DIR" "$TAG_DIR" "$PEAK_DIR" "$MACS3_DIR" "$HOMER_PEAK_DIR" "$PAUSING_DIR" "$LOG_DIR" "$MULTIQC_OUT"

check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
}

check_cmd_string() {
  local cmdstr="$1"
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
