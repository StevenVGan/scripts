#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# SEACR PEAK CALLING SCRIPT
#
# Runs 4 peak-calling passes from the same *_sorted.bam set:
#
# (A) WITH IgG CONTROL (IgG excluded as IP; first *IgG* BAM is control)
#   1) SEACR stringent -> ./peaks/seacr/stringent
#   2) SEACR relaxed   -> ./peaks/seacr/relaxed
#
# (B) NO CONTROL (includes IgG as IP too)
#   3) SEACR stringent -> ./peaks/seacr_2/stringent
#   4) SEACR relaxed   -> ./peaks/seacr_2/relaxed
#
# Each output directory gets a summary.tsv with:
#   sample, peak counts (raw + blacklist-filtered), FRiP (on filtered peaks)
#
# Workflow:
# 1) Clean BAM files (remove alt/random contigs) -> temporary clean BAMs
# 2) Generate bedGraph files from clean BAMs
# 3) Run SEACR with/without control, stringent/relaxed
# 4) Calculate FRiP statistics
# 5) Cleanup: remove clean BAMs, compress bedGraphs
###############################################################################

############################ CONFIG ###########################################

BASE="${HOME}/work/seq/CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka"

BAM_DIR="${BASE}/align/bam"
PEAK_DIR="${BASE}/peaks"

# Temporary directories for intermediate files
CLEAN_BAM_DIR="${BAM_DIR}/cleaned_bam"
BEDGRAPH_DIR="${BASE}/align/bedGraph"

# WITH IgG control outputs
SEACR_OUT_CTRL="${PEAK_DIR}/seacr"
SEACR_STRINGENT_CTRL="${SEACR_OUT_CTRL}/stringent"
SEACR_RELAXED_CTRL="${SEACR_OUT_CTRL}/relaxed"

# NO-control outputs
SEACR_OUT_NOCTRL="${PEAK_DIR}/seacr_2"
SEACR_STRINGENT_NOCTRL="${SEACR_OUT_NOCTRL}/stringent"
SEACR_RELAXED_NOCTRL="${SEACR_OUT_NOCTRL}/relaxed"

GENOME="hg38"
CHROM_SIZES="/mnt/share/archive/bkup/ref/genome/hg38/hg38.chrom.sizes"
BLACKLIST="${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed"

# SEACR settings
SEACR_THRESHOLD="0.01"    # Numeric threshold for peak calling
SEACR_NORM="non"          # For no-control: "non" or "norm"
SEACR_NORM_CTRL="norm"    # For with-control: "norm" or "non"

######################### END CONFIG ##########################################

# --- Logging (tee everything to ./logs) ---
LOG_DIR="${BASE}/logs"
mkdir -p "$LOG_DIR"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/3_peakcall_seacr_woCtrl_${RUN_TS}.log"
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
echo "For SEACR: SEACR_1.3.sh, samtools, bedtools, awk, gzip"
echo "========================="

check_cmd samtools
check_cmd bedtools
check_cmd awk
check_cmd gzip
check_cmd SEACR_1.3.sh

# Check for chrom.sizes file
if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "ERROR: Chromosome sizes file not found: $CHROM_SIZES" >&2
  exit 1
fi

mkdir -p "$CLEAN_BAM_DIR" "$BEDGRAPH_DIR"
mkdir -p "$SEACR_STRINGENT_CTRL" "$SEACR_RELAXED_CTRL"
mkdir -p "$SEACR_STRINGENT_NOCTRL" "$SEACR_RELAXED_NOCTRL"

count_lines() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo "0"
    return
  fi
  awk 'NF>0 && $0 !~ /^#/ {n++} END{print n+0}' "$f"
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

echo "=== SEACR peak calling run ==="
echo "[INFO] BASE:         $BASE"
echo "[INFO] BAM_DIR:      $BAM_DIR"
echo "[INFO] CLEAN_BAM:    $CLEAN_BAM_DIR"
echo "[INFO] BEDGRAPH:     $BEDGRAPH_DIR"
echo "[INFO] CHROM_SIZES:  $CHROM_SIZES"
echo "[INFO] BLACKLIST:    $BLACKLIST"
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
IGG_CONTROL_BEDGRAPH=""
if (( ${#igg_bams[@]} > 0 )); then
  IGG_CONTROL_BAM="${igg_bams[0]}"
  IGG_CONTROL_SAMPLE="$(basename "$IGG_CONTROL_BAM" _sorted.bam)"
  IGG_CONTROL_BEDGRAPH="${BEDGRAPH_DIR}/${IGG_CONTROL_SAMPLE}_cleaned.bedgraph"
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
# STEP 1: Clean BAMs and Generate bedGraphs for ALL samples
###############################################################################
echo "=== STEP 1: Cleaning BAMs and generating bedGraphs ==="

for bam in "${bams[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"
  
  clean_bam="${CLEAN_BAM_DIR}/${sample}_cleaned.bam"
  bedgraph="${BEDGRAPH_DIR}/${sample}_cleaned.bedgraph"
  
  echo "[Processing] $sample"
  
  # Clean BAM: remove alt/random contigs
  samtools view -h "$bam" \
    | awk '!/chr.*_alt|chr.*_random/' \
    | samtools view -b - > "$clean_bam"
  
  # Generate bedGraph
  bedtools bamtobed -i "$clean_bam" \
    | bedtools genomecov -bg -i - -g "$CHROM_SIZES" \
    > "$bedgraph"
  
  echo "  -> Clean BAM:  $clean_bam"
  echo "  -> bedGraph:   $bedgraph"
done

echo "[DONE] All bedGraphs generated"
echo

###############################################################################
# (A1) SEACR WITH IgG CONTROL - STRINGENT
###############################################################################
SEACR_STRINGENT_CTRL_TSV="${SEACR_STRINGENT_CTRL}/summary.tsv"
echo -e "sample\tmode\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$SEACR_STRINGENT_CTRL_TSV"

if [[ -n "$IGG_CONTROL_BAM" ]]; then
  echo "=== (A1) SEACR WITH IgG CONTROL - STRINGENT ==="
  
  for bam in "${ip_bams_ctrl[@]}"; do
    base=$(basename "$bam")
    sample="${base%_sorted.bam}"
    bedgraph="${BEDGRAPH_DIR}/${sample}_cleaned.bedgraph"
    
    echo "[SEACR stringent ctrl] $sample  (control: $IGG_CONTROL_SAMPLE)"
    
    output_prefix="${SEACR_STRINGENT_CTRL}/${sample}_vsIgG"
    
    SEACR_1.3.sh \
      "$bedgraph" \
      "$IGG_CONTROL_BEDGRAPH" \
      "$SEACR_NORM_CTRL" \
      stringent \
      "$output_prefix"
    
    # SEACR outputs: {prefix}.stringent.bed
    peak_file="${output_prefix}.stringent.bed"
    
    if [[ ! -f "$peak_file" ]]; then
      echo -e "${sample}\tstringent_ctrl\t0\t0\t0\t0\t0.0000" >> "$SEACR_STRINGENT_CTRL_TSV"
      continue
    fi
    
    n_raw=$(count_lines "$peak_file")
    
    filtered_bed="${output_prefix}_filtered.bed"
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
    
    echo -e "${sample}\tstringent_ctrl\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$SEACR_STRINGENT_CTRL_TSV"
  done
  
  echo "[DONE] SEACR stringent WITH-control summary: $SEACR_STRINGENT_CTRL_TSV"
  echo
fi

###############################################################################
# (A2) SEACR WITH IgG CONTROL - RELAXED
###############################################################################
SEACR_RELAXED_CTRL_TSV="${SEACR_RELAXED_CTRL}/summary.tsv"
echo -e "sample\tmode\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$SEACR_RELAXED_CTRL_TSV"

if [[ -n "$IGG_CONTROL_BAM" ]]; then
  echo "=== (A2) SEACR WITH IgG CONTROL - RELAXED ==="
  
  for bam in "${ip_bams_ctrl[@]}"; do
    base=$(basename "$bam")
    sample="${base%_sorted.bam}"
    bedgraph="${BEDGRAPH_DIR}/${sample}_cleaned.bedgraph"
    
    echo "[SEACR relaxed ctrl] $sample  (control: $IGG_CONTROL_SAMPLE)"
    
    output_prefix="${SEACR_RELAXED_CTRL}/${sample}_vsIgG"
    
    SEACR_1.3.sh \
      "$bedgraph" \
      "$IGG_CONTROL_BEDGRAPH" \
      "$SEACR_NORM_CTRL" \
      relaxed \
      "$output_prefix"
    
    # SEACR outputs: {prefix}.relaxed.bed
    peak_file="${output_prefix}.relaxed.bed"
    
    if [[ ! -f "$peak_file" ]]; then
      echo -e "${sample}\trelaxed_ctrl\t0\t0\t0\t0\t0.0000" >> "$SEACR_RELAXED_CTRL_TSV"
      continue
    fi
    
    n_raw=$(count_lines "$peak_file")
    
    filtered_bed="${output_prefix}_filtered.bed"
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
    
    echo -e "${sample}\trelaxed_ctrl\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$SEACR_RELAXED_CTRL_TSV"
  done
  
  echo "[DONE] SEACR relaxed WITH-control summary: $SEACR_RELAXED_CTRL_TSV"
  echo
fi

###############################################################################
# (B1) SEACR NO CONTROL - STRINGENT
###############################################################################
SEACR_STRINGENT_NOCTRL_TSV="${SEACR_STRINGENT_NOCTRL}/summary.tsv"
echo -e "sample\tmode\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$SEACR_STRINGENT_NOCTRL_TSV"

echo "=== (B1) SEACR NO CONTROL - STRINGENT (includes IgG) ==="

for bam in "${ip_bams_noctrl[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"
  bedgraph="${BEDGRAPH_DIR}/${sample}_cleaned.bedgraph"
  
  echo "[SEACR stringent noc] $sample"
  
  output_prefix="${SEACR_STRINGENT_NOCTRL}/${sample}_noctrl"
  
  SEACR_1.3.sh \
    "$bedgraph" \
    "$SEACR_THRESHOLD" \
    "$SEACR_NORM" \
    stringent \
    "$output_prefix"
  
  # SEACR outputs: {prefix}.stringent.bed
  peak_file="${output_prefix}.stringent.bed"
  
  if [[ ! -f "$peak_file" ]]; then
    echo -e "${sample}\tstringent_noctrl\t0\t0\t0\t0\t0.0000" >> "$SEACR_STRINGENT_NOCTRL_TSV"
    continue
  fi
  
  n_raw=$(count_lines "$peak_file")
  
  filtered_bed="${output_prefix}_filtered.bed"
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
  
  echo -e "${sample}\tstringent_noctrl\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$SEACR_STRINGENT_NOCTRL_TSV"
done

echo "[DONE] SEACR stringent NO-control summary: $SEACR_STRINGENT_NOCTRL_TSV"
echo

###############################################################################
# (B2) SEACR NO CONTROL - RELAXED
###############################################################################
SEACR_RELAXED_NOCTRL_TSV="${SEACR_RELAXED_NOCTRL}/summary.tsv"
echo -e "sample\tmode\tn_peaks_raw\tn_peaks_filtered\treads_in_peaks\ttotal_mapped\tfrip" > "$SEACR_RELAXED_NOCTRL_TSV"

echo "=== (B2) SEACR NO CONTROL - RELAXED (includes IgG) ==="

for bam in "${ip_bams_noctrl[@]}"; do
  base=$(basename "$bam")
  sample="${base%_sorted.bam}"
  bedgraph="${BEDGRAPH_DIR}/${sample}_cleaned.bedgraph"
  
  echo "[SEACR relaxed noc] $sample"
  
  output_prefix="${SEACR_RELAXED_NOCTRL}/${sample}_noctrl"
  
  SEACR_1.3.sh \
    "$bedgraph" \
    "$SEACR_THRESHOLD" \
    "$SEACR_NORM" \
    relaxed \
    "$output_prefix"
  
  # SEACR outputs: {prefix}.relaxed.bed
  peak_file="${output_prefix}.relaxed.bed"
  
  if [[ ! -f "$peak_file" ]]; then
    echo -e "${sample}\trelaxed_noctrl\t0\t0\t0\t0\t0.0000" >> "$SEACR_RELAXED_NOCTRL_TSV"
    continue
  fi
  
  n_raw=$(count_lines "$peak_file")
  
  filtered_bed="${output_prefix}_filtered.bed"
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
  
  echo -e "${sample}\trelaxed_noctrl\t${n_raw}\t${n_filt}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$SEACR_RELAXED_NOCTRL_TSV"
done

echo "[DONE] SEACR relaxed NO-control summary: $SEACR_RELAXED_NOCTRL_TSV"
echo

###############################################################################
# CLEANUP: Remove clean BAMs and compress bedGraphs
###############################################################################
echo "=== CLEANUP: Removing clean BAMs and compressing bedGraphs ==="

# Remove clean BAMs
echo "[INFO] Removing clean BAMs from $CLEAN_BAM_DIR"
rm -f "$CLEAN_BAM_DIR"/*.bam
rmdir "$CLEAN_BAM_DIR" 2>/dev/null || true

# Compress bedGraphs
echo "[INFO] Compressing bedGraphs in $BEDGRAPH_DIR"
for bg in "$BEDGRAPH_DIR"/*.bedgraph; do
  if [[ -f "$bg" ]]; then
    echo "  Compressing: $(basename "$bg")"
    gzip -f "$bg"
  fi
done

echo "[DONE] Cleanup complete"
echo

echo "=== All done ==="
echo "WITH-control stringent: $SEACR_STRINGENT_CTRL (summary.tsv)"
echo "WITH-control relaxed:   $SEACR_RELAXED_CTRL (summary.tsv)"
echo "NO-control stringent:   $SEACR_STRINGENT_NOCTRL (summary.tsv)"
echo "NO-control relaxed:     $SEACR_RELAXED_NOCTRL (summary.tsv)"
echo
echo "Clean BAMs removed from: $CLEAN_BAM_DIR"
echo "bedGraphs compressed in: $BEDGRAPH_DIR"