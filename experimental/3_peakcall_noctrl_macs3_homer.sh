#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# Peak calling "round 2" (NO CONTROL) on all BAMs:
#   - MACS3 callpeak for every *_sorted.bam (including IgG)
#   - HOMER findPeaks for every *_sorted.bam (including IgG), reusing existing
#     tag directories from align/tags when present (fallback to build if missing)
#
# Outputs (NEW dirs only; does NOT modify prior results):
#   peaks/macs3_2/   + summary.tsv (FRiP + peak counts)
#   peaks/homer_2/   + summary.tsv (FRiP + peak counts)
###############################################################################

############################ CONFIG ###########################################

# Project base directory (match your pipeline BASE)
BASE="${HOME}/work/seq/CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka"

# Input BAMs and existing HOMER tagdirs (from upstream pipeline)
BAM_DIR="${BASE}/align/bam"
TAG_DIR_EXISTING="${BASE}/align/tags"

# Output dirs (NEW)
PEAK_DIR="${BASE}/peaks"
MACS3_OUT="${PEAK_DIR}/macs3_2"
HOMER_OUT="${PEAK_DIR}/homer_2"

# Genome / blacklist (optional but recommended; keep consistent with your pipeline)
GENOME="hg38"
BLACKLIST="${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed"

# MACS3 parameters (no control)
MACS3_GENOMESIZE="hs"
MACS3_FORMAT="BAMPE"     # change to BAM if single-end
MACS3_Q=0.05
MACS3_KEEP_DUP="auto"    # try "all" if peak counts are suspiciously low

# How to run macs3:
#   - If macs3 is on PATH: set USE_CONDA_MACS3=0
#   - If you still rely on a conda env named macs3: set USE_CONDA_MACS3=1
USE_CONDA_MACS3=0

# HOMER peak calling style:
#   "factor" for TF-like narrow peaks; "histone" for broad marks
HOMER_STYLE="factor"

######################### END CONFIG ##########################################

check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
}

# Requirements
check_cmd samtools
check_cmd bedtools
check_cmd findPeaks
check_cmd pos2bed.pl
if [[ "$USE_CONDA_MACS3" -eq 1 ]]; then
  check_cmd conda
else
  check_cmd macs3
fi
# makeTagDirectory only needed for fallback when tagdir is missing
if command -v makeTagDirectory >/dev/null 2>&1; then
  HAVE_MAKETAG=1
else
  HAVE_MAKETAG=0
fi

mkdir -p "$MACS3_OUT" "$HOMER_OUT"

# Choose macs3 invocation
if [[ "$USE_CONDA_MACS3" -eq 1 ]]; then
  MACS3_CMD=( conda run -n macs3 macs3 )
else
  MACS3_CMD=( macs3 )
fi

# Count non-empty, non-comment lines
count_lines() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo "0"
    return
  fi
  awk 'NF>0 && $0 !~ /^#/ {n++} END{print n+0}' "$f"
}

# HOMER peaks.txt typically has a header line starting with "PeakID"
count_homer_peaks_txt() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo "0"
    return
  fi
  awk '
    BEGIN{n=0; started=0}
    /^#/ {next}
    {
      if (started==0) {
        started=1
        if ($1=="PeakID") next
      }
      if (NF>0) n++
    }
    END{print n+0}
  ' "$f"
}

# FRiP calculator: BAM + BED -> reads_in_peaks, total_mapped, frip
# total_mapped excludes unmapped(0x4) + secondary(0x100) + supplementary(0x800) = 2308
calc_frip() {
  local bam="$1"
  local bed="$2"

  local total_mapped reads_in_peaks frip
  total_mapped=$(samtools view -c -F 2308 "$bam" || echo 0)

  if [[ ! -f "$bed" ]]; then
    reads_in_peaks=0
    frip="0.0000"
    echo -e "${reads_in_peaks}\t${total_mapped}\t${frip}"
    return
  fi

  reads_in_peaks=$(
    bedtools intersect -a "$bam" -b "$bed" -u \
      | samtools view -c - \
      || echo 0
  )

  frip=$(awk -v inpeaks="$reads_in_peaks" -v total="$total_mapped" \
    'BEGIN { if (total > 0) printf "%.4f", inpeaks/total; else print "0.0000" }')

  echo -e "${reads_in_peaks}\t${total_mapped}\t${frip}"
}

echo "=== Peak calling round 2 (NO CONTROL) ==="
echo "[INFO] BASE:            $BASE"
echo "[INFO] BAM_DIR:         $BAM_DIR"
echo "[INFO] TAG_DIR_EXISTING $TAG_DIR_EXISTING"
echo "[INFO] MACS3_OUT:       $MACS3_OUT"
echo "[INFO] HOMER_OUT:       $HOMER_OUT"
echo "[INFO] BLACKLIST:       $BLACKLIST"

bams=( "$BAM_DIR"/*_sorted.bam )
if (( ${#bams[@]} == 0 )); then
  echo "ERROR: No *_sorted.bam found in $BAM_DIR"
  exit 1
fi

###############################################################################
# 1) MACS3 NO-CONTROL on all BAMs
###############################################################################
MACS3_TSV="${MACS3_OUT}/summary.tsv"
echo -e "sample\tpeak_type\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip\tpeak_file\tfiltered_bed" > "$MACS3_TSV"

for bam in "${bams[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"
  name="${sample}_noctrl"

  echo "=== [MACS3] $sample ==="

  "${MACS3_CMD[@]}" callpeak \
    -t "$bam" \
    -g "$MACS3_GENOMESIZE" \
    -q "$MACS3_Q" \
    --format "$MACS3_FORMAT" \
    --keep-dup "$MACS3_KEEP_DUP" \
    --outdir "$MACS3_OUT" \
    --name "$name"

  narrow="${MACS3_OUT}/${name}_peaks.narrowPeak"
  broad="${MACS3_OUT}/${name}_peaks.broadPeak"

  peak_type=""
  peak_file=""
  if [[ -f "$narrow" ]]; then
    peak_type="narrowPeak"
    peak_file="$narrow"
  elif [[ -f "$broad" ]]; then
    peak_type="broadPeak"
    peak_file="$broad"
  else
    echo "[WARN] MACS3 produced no narrow/broad peak file for $sample"
    echo -e "${sample}\tNA\t0\t0\t0\t0\t0.0000\tNA\tNA" >> "$MACS3_TSV"
    continue
  fi

  n_raw=$(count_lines "$peak_file")

  filtered_bed="${MACS3_OUT}/${name}_filtered.bed"
  if [[ -f "$BLACKLIST" ]]; then
    bedtools intersect -a "$peak_file" -b "$BLACKLIST" -v > "$filtered_bed"
  else
    cp -f "$peak_file" "$filtered_bed"
  fi
  n_filt=$(count_lines "$filtered_bed")

  frip_vals=$(calc_frip "$bam" "$filtered_bed")
  reads_in_peaks=$(echo "$frip_vals" | cut -f1)
  total_mapped=$(echo "$frip_vals" | cut -f2)
  frip=$(echo "$frip_vals" | cut -f3)

  echo -e "${sample}\t${peak_type}\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}\t${peak_file}\t${filtered_bed}" >> "$MACS3_TSV"
done

echo "[DONE] MACS3 summary: $MACS3_TSV"

###############################################################################
# 2) HOMER NO-CONTROL on all BAMs (reuse align/tags when present)
###############################################################################
HOMER_TSV="${HOMER_OUT}/summary.tsv"
echo -e "sample\thomer_style\tn_peaks_txt\tn_peaks_bed_filtered\treads_in_peaks\ttotal_mapped\tfrip\tpeaks_txt\tpeaks_bed_filtered\ttagdir_used" > "$HOMER_TSV"

for bam in "${bams[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"

  echo "=== [HOMER] $sample ==="

  # Prefer existing tagdir from upstream pipeline
  tagdir="${TAG_DIR_EXISTING}/${sample}"

  # Fallback: build a tagdir under homer_2 only if missing
  if [[ ! -d "$tagdir" ]]; then
    echo "[WARN] Tagdir not found at $tagdir; attempting to create from BAM"
    if [[ "$HAVE_MAKETAG" -ne 1 ]]; then
      echo "ERROR: makeTagDirectory not found, but tagdir missing for $sample. Install HOMER or add it to PATH." >&2
      exit 1
    fi
    mkdir -p "${HOMER_OUT}/tagdirs"
    tagdir="${HOMER_OUT}/tagdirs/${sample}"
    makeTagDirectory "$tagdir" "$bam" -format bam -tbp 1
  fi

  peaks_txt="${HOMER_OUT}/${sample}.peaks.txt"
  peaks_bed="${HOMER_OUT}/${sample}.peaks.bed"
  peaks_bed_filt="${HOMER_OUT}/${sample}.peaks.filtered.bed"

  findPeaks "$tagdir" -style "$HOMER_STYLE" -o "$peaks_txt"
  n_txt=$(count_homer_peaks_txt "$peaks_txt")

  # Convert HOMER output to BED for FRiP
  pos2bed.pl "$peaks_txt" > "$peaks_bed"

  # Blacklist filter
  if [[ -f "$BLACKLIST" ]]; then
    bedtools intersect -a "$peaks_bed" -b "$BLACKLIST" -v > "$peaks_bed_filt"
  else
    cp -f "$peaks_bed" "$peaks_bed_filt"
  fi
  n_bed_filt=$(count_lines "$peaks_bed_filt")

  frip_vals=$(calc_frip "$bam" "$peaks_bed_filt")
  reads_in_peaks=$(echo "$frip_vals" | cut -f1)
  total_mapped=$(echo "$frip_vals" | cut -f2)
  frip=$(echo "$frip_vals" | cut -f3)

  echo -e "${sample}\t${HOMER_STYLE}\t${n_txt}\t${n_bed_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}\t${peaks_txt}\t${peaks_bed_filt}\t${tagdir}" >> "$HOMER_TSV"
done

echo "[DONE] HOMER summary: $HOMER_TSV"
echo "=== Finished peak calling round 2 (NO CONTROL) ==="
