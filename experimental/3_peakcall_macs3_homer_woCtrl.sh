#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# INSPECTIVE PEAK CALLING SCRIPT
#
# Runs 4 peak-calling passes from the same *_sorted.bam set:
#
# (A) WITH IgG CONTROL (IgG excluded as IP; first *IgG* BAM is control)
#   1) MACS3 callpeak  -> ./peaks/macs3
#   2) HOMER findPeaks -> ./peaks/homer   (uses IgG tagdir as -i background)
#
# (B) NO CONTROL (includes IgG as IP too)
#   3) MACS3 callpeak  -> ./peaks/macs3_2
#   4) HOMER findPeaks -> ./peaks/homer_2
#
# Each output directory gets a summary.tsv with:
#   sample, peak counts (raw + blacklist-filtered), FRiP (on filtered peaks)
#
# Reuses existing HOMER tag directories at: align/tags/<sample>
# Falls back to makeTagDirectory only if a required tagdir is missing.
###############################################################################

############################ CONFIG ###########################################

BASE="${HOME}/work/seq/CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka"

BAM_DIR="${BASE}/align/bam"
TAG_DIR_EXISTING="${BASE}/align/tags"

PEAK_DIR="${BASE}/peaks"

# WITH IgG control outputs
MACS3_OUT_CTRL="${PEAK_DIR}/macs3"
HOMER_OUT_CTRL="${PEAK_DIR}/homer"

# NO-control outputs
MACS3_OUT_NOCTRL="${PEAK_DIR}/macs3_2"
HOMER_OUT_NOCTRL="${PEAK_DIR}/homer_2"

GENOME="hg38"
BLACKLIST="${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed"

# MACS3 settings
MACS3_GENOMESIZE="hs"
MACS3_FORMAT="BAMPE"     # change to BAM if single-end
MACS3_Q=0.05
MACS3_KEEP_DUP="auto"    # try "all" if you suspect duplicates handling is killing peaks

# How to run macs3
USE_CONDA_MACS3=0

# HOMER style
HOMER_STYLE="factor"     # or "histone"

######################### END CONFIG ##########################################

# --- Logging (tee everything to ./logs) ---
LOG_DIR="${BASE}/logs"
mkdir -p "$LOG_DIR"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/3_peakcall_macs3_homer_woCtrl_${RUN_TS}.log"
echo "[INFO] Logging to: $LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1
# --- end logging ---

check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
}

echo "=== Required programs ==="
echo "General: bash, coreutils"
echo "For MACS3: macs3 (or conda), samtools, bedtools"
echo "For HOMER: findPeaks, pos2bed.pl (and makeTagDirectory only if tagdirs missing)"
echo "========================="

check_cmd samtools
check_cmd bedtools
check_cmd findPeaks
check_cmd pos2bed.pl
if [[ "$USE_CONDA_MACS3" -eq 1 ]]; then
  check_cmd conda
else
  check_cmd macs3
fi

HAVE_MAKETAG=0
if command -v makeTagDirectory >/dev/null 2>&1; then
  HAVE_MAKETAG=1
fi

mkdir -p "$MACS3_OUT_CTRL" "$HOMER_OUT_CTRL" "$MACS3_OUT_NOCTRL" "$HOMER_OUT_NOCTRL"

if [[ "$USE_CONDA_MACS3" -eq 1 ]]; then
  MACS3_CMD=( conda run -n macs3 macs3 )
else
  MACS3_CMD=( macs3 )
fi

count_lines() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo "0"
    return
  fi
  awk 'NF>0 && $0 !~ /^#/ {n++} END{print n+0}' "$f"
}

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

# FRiP: total_mapped excludes unmapped(0x4) + secondary(0x100) + supplementary(0x800) = 2308
calc_frip() {
  local bam="$1"
  local bed="$2"

  local total_mapped reads_in_peaks frip
  total_mapped=$(samtools view -c -F 2308 "$bam" || echo 0)

  if [[ ! -f "$bed" ]]; then
    echo -e "0\t${total_mapped}\t0.0000"
    return
  fi

  reads_in_peaks=$(
    bedtools intersect -abam "$bam" -b "$bed" -u \
      | samtools view -c - \
      || echo 0
  )

  frip=$(awk -v inpeaks="$reads_in_peaks" -v total="$total_mapped" \
    'BEGIN { if (total > 0) printf "%.4f", inpeaks/total; else print "0.0000" }')

  echo -e "${reads_in_peaks}\t${total_mapped}\t${frip}"
}

resolve_tagdir_or_build() {
  # args: sample bam out_tagdir_parent
  local sample="$1"
  local bam="$2"
  local out_parent="$3"

  local tagdir="${TAG_DIR_EXISTING}/${sample}"
  if [[ -d "$tagdir" ]]; then
    echo "$tagdir"
    return
  fi

  if [[ "$HAVE_MAKETAG" -ne 1 ]]; then
    echo "ERROR: Missing HOMER tagdir for sample '$sample' at '$tagdir' and makeTagDirectory not found." >&2
    exit 1
  fi

  mkdir -p "$out_parent"
  tagdir="${out_parent}/${sample}"
  makeTagDirectory "$tagdir" "$bam" -format bam -tbp 1
  echo "$tagdir"
}

echo "=== Peak calling inspect run ==="
echo "[INFO] BASE:            $BASE"
echo "[INFO] BAM_DIR:         $BAM_DIR"
echo "[INFO] TAG_DIR_EXISTING $TAG_DIR_EXISTING"
echo "[INFO] BLACKLIST:       $BLACKLIST"
echo

bams=( "$BAM_DIR"/*_sorted.bam )
if (( ${#bams[@]} == 0 )); then
  echo "ERROR: No *_sorted.bam found in $BAM_DIR"
  exit 1
fi

# Identify IgG control BAM (simple heuristic: first BAM filename containing "IgG" case-insensitive)
igg_bams=()
for b in "${bams[@]}"; do
  bn=$(basename "$b")
  if [[ "$bn" =~ [Ii][Gg][Gg] ]]; then
    igg_bams+=( "$b" )
  fi
done

IGG_CONTROL_BAM=""
IGG_CONTROL_SAMPLE=""
if (( ${#igg_bams[@]} > 0 )); then
  IGG_CONTROL_BAM="${igg_bams[0]}"
  IGG_CONTROL_SAMPLE="$(basename "$IGG_CONTROL_BAM" _sorted.bam)"
  if (( ${#igg_bams[@]} > 1 )); then
    echo "[WARN] Multiple IgG BAMs detected; using FIRST as control:"
    echo "       $IGG_CONTROL_BAM"
  else
    echo "[INFO] IgG control BAM:"
    echo "       $IGG_CONTROL_BAM"
  fi
else
  echo "[WARN] No IgG BAM detected (filename contains 'IgG'); WITH-control runs will be skipped."
fi
echo

# Build IP lists:
ip_bams_noctrl=( "${bams[@]}" )

ip_bams_ctrl=()
for b in "${bams[@]}"; do
  bn=$(basename "$b")
  if [[ "$bn" =~ [Ii][Gg][Gg] ]]; then
    continue  # exclude IgG as IP for WITH-control runs
  fi
  ip_bams_ctrl+=( "$b" )
done

###############################################################################
# (A1) MACS3 WITH IgG CONTROL  -> peaks/macs3
###############################################################################
MACS3_CTRL_TSV="${MACS3_OUT_CTRL}/summary.tsv"
echo -e "sample\tpeak_type\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$MACS3_CTRL_TSV"

if [[ -n "$IGG_CONTROL_BAM" ]]; then
  echo "=== (A1) MACS3 WITH IgG CONTROL ==="
  for bam in "${ip_bams_ctrl[@]}"; do
    base=$(basename "$bam")
    sample="${base%_sorted.bam}"
    name="${sample}_vsIgG"

    echo "[MACS3 ctrl] $sample  (control: $IGG_CONTROL_SAMPLE)"

    "${MACS3_CMD[@]}" callpeak \
      -t "$bam" \
      -c "$IGG_CONTROL_BAM" \
      -g "$MACS3_GENOMESIZE" \
      -q "$MACS3_Q" \
      --format "$MACS3_FORMAT" \
      --keep-dup "$MACS3_KEEP_DUP" \
      --outdir "$MACS3_OUT_CTRL" \
      --name "$name"

    narrow="${MACS3_OUT_CTRL}/${name}_peaks.narrowPeak"
    broad="${MACS3_OUT_CTRL}/${name}_peaks.broadPeak"

    peak_type="NA"
    peak_file=""
    if [[ -f "$narrow" ]]; then
      peak_type="narrowPeak"
      peak_file="$narrow"
    elif [[ -f "$broad" ]]; then
      peak_type="broadPeak"
      peak_file="$broad"
    fi

    if [[ -z "$peak_file" ]]; then
      echo -e "${sample}\tNA\t0\t0\t0\t0\t0.0000" >> "$MACS3_CTRL_TSV"
      continue
    fi

    n_raw=$(count_lines "$peak_file")

    filtered_bed="${MACS3_OUT_CTRL}/${name}_filtered.bed"
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

    echo -e "${sample}\t${peak_type}\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$MACS3_CTRL_TSV"
  done
  echo "[DONE] MACS3 WITH-control summary: $MACS3_CTRL_TSV"
  echo
fi

###############################################################################
# (A2) HOMER WITH IgG BACKGROUND -> peaks/homer
###############################################################################
HOMER_CTRL_TSV="${HOMER_OUT_CTRL}/summary.tsv"
echo -e "sample\thomer_style\tn_peaks_txt\tn_peaks_bed_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$HOMER_CTRL_TSV"

if [[ -n "$IGG_CONTROL_BAM" ]]; then
  echo "=== (A2) HOMER WITH IgG BACKGROUND ==="

  # Resolve/control tagdir (reuse align/tags; fallback build under peaks/homer/_control_tagdir)
  control_tagdir="$(resolve_tagdir_or_build "$IGG_CONTROL_SAMPLE" "$IGG_CONTROL_BAM" "${HOMER_OUT_CTRL}/_control_tagdir")"
  echo "[INFO] HOMER control tagdir: $control_tagdir"
  echo

  for bam in "${ip_bams_ctrl[@]}"; do
    base=$(basename "$bam")
    sample="${base%_sorted.bam}"

    echo "[HOMER ctrl] $sample  (background: $IGG_CONTROL_SAMPLE)"

    tagdir="$(resolve_tagdir_or_build "$sample" "$bam" "${HOMER_OUT_CTRL}/_tagdirs_fallback")"

    peaks_txt="${HOMER_OUT_CTRL}/${sample}_vsIgG.peaks.txt"
    peaks_bed="${HOMER_OUT_CTRL}/${sample}_vsIgG.peaks.bed"
    peaks_bed_filt="${HOMER_OUT_CTRL}/${sample}_vsIgG.peaks.filtered.bed"

    findPeaks "$tagdir" -style "$HOMER_STYLE" -i "$control_tagdir" -o "$peaks_txt"
    n_txt=$(count_homer_peaks_txt "$peaks_txt")

    pos2bed.pl "$peaks_txt" > "$peaks_bed"

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

    echo -e "${sample}\t${HOMER_STYLE}\t${n_txt}\t${n_bed_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$HOMER_CTRL_TSV"
  done

  echo "[DONE] HOMER WITH-control summary: $HOMER_CTRL_TSV"
  echo
fi

###############################################################################
# (B1) MACS3 NO CONTROL -> peaks/macs3_2  (includes IgG as IP)
###############################################################################
MACS3_NOCTRL_TSV="${MACS3_OUT_NOCTRL}/summary.tsv"
echo -e "sample\tpeak_type\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$MACS3_NOCTRL_TSV"

echo "=== (B1) MACS3 NO CONTROL (includes IgG) ==="
for bam in "${ip_bams_noctrl[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"
  name="${sample}_noctrl"

  echo "[MACS3 noc] $sample"

  "${MACS3_CMD[@]}" callpeak \
    -t "$bam" \
    -g "$MACS3_GENOMESIZE" \
    -q "$MACS3_Q" \
    --format "$MACS3_FORMAT" \
    --keep-dup "$MACS3_KEEP_DUP" \
    --outdir "$MACS3_OUT_NOCTRL" \
    --name "$name"

  narrow="${MACS3_OUT_NOCTRL}/${name}_peaks.narrowPeak"
  broad="${MACS3_OUT_NOCTRL}/${name}_peaks.broadPeak"

  peak_type="NA"
  peak_file=""
  if [[ -f "$narrow" ]]; then
    peak_type="narrowPeak"
    peak_file="$narrow"
  elif [[ -f "$broad" ]]; then
    peak_type="broadPeak"
    peak_file="$broad"
  fi

  if [[ -z "$peak_file" ]]; then
    echo -e "${sample}\tNA\t0\t0\t0\t0\t0.0000" >> "$MACS3_NOCTRL_TSV"
    continue
  fi

  n_raw=$(count_lines "$peak_file")

  filtered_bed="${MACS3_OUT_NOCTRL}/${name}_filtered.bed"
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

  echo -e "${sample}\t${peak_type}\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$MACS3_NOCTRL_TSV"
done
echo "[DONE] MACS3 NO-control summary: $MACS3_NOCTRL_TSV"
echo

###############################################################################
# (B2) HOMER NO CONTROL -> peaks/homer_2  (includes IgG as IP)
###############################################################################
HOMER_NOCTRL_TSV="${HOMER_OUT_NOCTRL}/summary.tsv"
echo -e "sample\thomer_style\tn_peaks_txt\tn_peaks_bed_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$HOMER_NOCTRL_TSV"

echo "=== (B2) HOMER NO CONTROL (includes IgG) ==="
for bam in "${ip_bams_noctrl[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"

  echo "[HOMER noc] $sample"

  tagdir="$(resolve_tagdir_or_build "$sample" "$bam" "${HOMER_OUT_NOCTRL}/_tagdirs_fallback")"

  peaks_txt="${HOMER_OUT_NOCTRL}/${sample}_noctrl.peaks.txt"
  peaks_bed="${HOMER_OUT_NOCTRL}/${sample}_noctrl.peaks.bed"
  peaks_bed_filt="${HOMER_OUT_NOCTRL}/${sample}_noctrl.peaks.filtered.bed"

  findPeaks "$tagdir" -style "$HOMER_STYLE" -o "$peaks_txt"
  n_txt=$(count_homer_peaks_txt "$peaks_txt")

  pos2bed.pl "$peaks_txt" > "$peaks_bed"

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

  echo -e "${sample}\t${HOMER_STYLE}\t${n_txt}\t${n_bed_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$HOMER_NOCTRL_TSV"
done
echo "[DONE] HOMER NO-control summary: $HOMER_NOCTRL_TSV"
echo

echo "=== All done ==="
echo "WITH-control:  $MACS3_OUT_CTRL (summary.tsv), $HOMER_OUT_CTRL (summary.tsv)"
echo "NO-control:    $MACS3_OUT_NOCTRL (summary.tsv), $HOMER_OUT_NOCTRL (summary.tsv)"
