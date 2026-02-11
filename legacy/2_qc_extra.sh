#!/usr/bin/env bash
set -euo pipefail

# TODO: TSS annotation is a mess. Fix later.
# TODO: The Multiqc is still a bit messy. 1) SPP Results are not in the section; 2) Strand shift is not plotted; 3) Order of general stats is weird, so does the section order; 4) HOMER chromatin annotation section could be improved.
# TODO: Mix this script with the upstream one?

###############################################################################
# CUT&RUN EXTRA QC + MULTIQC SCRIPT
#
# Run this AFTER the upstream pipeline.
###############################################################################

############################ CONFIG SECTION ###################################

BASE="${HOME}/work/seq/CUTRUN/test"
GENOME="hg38"
ANNOT_BED="${BASE}/ref/hg38_tss_regions.bed"

# Directories
ALIGN_DIR="${BASE}/align"
BAM_DIR="${ALIGN_DIR}/bam"
BAMQC_DIR="${ALIGN_DIR}/bamqc"
PEAK_DIR="${BASE}/peaks"
MACS2_DIR="${PEAK_DIR}/macs2"
HOMER_ANN_DIR="${PEAK_DIR}/homer"
LOG_DIR="${BASE}/logs"
MULTIQC_OUT="${BASE}/multiqc"

# Pre-calculated FRiP table from upstream script
FRIP_TSV="${MACS2_DIR}/frip_scores.tsv"
THREADS=8

######################## END CONFIG SECTION ###################################

check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
}

echo "=== CUT&RUN extra QC + MultiQC starting ==="

for c in samtools preseq run_spp.R plotFingerprint computeMatrix plotProfile \
         plotPCA plotCorrelation annotatePeaks.pl multiqc; do
  check_cmd "$c"
done

mkdir -p "$BAMQC_DIR" "$HOMER_ANN_DIR" "$MULTIQC_OUT" "$LOG_DIR"

QC_LOG="${LOG_DIR}/2_cutrun_qc_extra_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to: $QC_LOG"
exec > >(tee -a "$QC_LOG") 2>&1

###############################################################################
# STEP 1: Per-BAM QC
###############################################################################

echo "=== STEP 1: Per-BAM QC ==="
sorted_bams=( "${BAM_DIR}"/*_sorted.bam )
if (( ${#sorted_bams[@]} == 0 )); then
  echo "[WARN] No *_sorted.bam files found in $BAM_DIR. Exiting."
  exit 0
fi

for bam in "${sorted_bams[@]}"; do
  sample=$(basename "$bam" _sorted.bam)
  echo "[STEP1] Sample: $sample"

  # 1.1 samtools stats
  stats_file="${BAMQC_DIR}/${sample}.stats"
  if [[ ! -f "$stats_file" ]]; then
    echo "  - samtools stats -> $stats_file"
    samtools stats "$bam" > "$stats_file"
  fi

  # 1.2 preseq
  preseq_file="${BAMQC_DIR}/${sample}.preseq.txt"
  if [[ ! -f "$preseq_file" ]]; then
    echo "  - preseq lc_extrap -> $preseq_file"
    preseq lc_extrap -B "$bam" -o "$preseq_file" -v
  fi

  # 1.3 phantompeakqualtools
  spp_out="${BAMQC_DIR}/${sample}.spp.out"
  spp_pdf="${BAMQC_DIR}/${sample}.spp.pdf"
  if [[ ! -f "$spp_out" ]]; then
    echo "  - run_spp.R -> $spp_out"
    run_spp.R -c="$bam" -p="$THREADS" -savp="$spp_pdf" -out="$spp_out"
  fi
done

###############################################################################
# STEP 2: deepTools QC
###############################################################################

echo "=== STEP 2: deepTools QC ==="
FINGERPRINT_RAW="${BAMQC_DIR}/fingerprint_rawcounts.txt"
FINGERPRINT_PNG="${BAMQC_DIR}/fingerprint.png"

if [[ ! -f "$FINGERPRINT_RAW" ]]; then
  plotFingerprint -b "${sorted_bams[@]}" -o "$FINGERPRINT_PNG" \
    --outRawCounts "$FINGERPRINT_RAW" --ignoreDuplicates -p "$THREADS"
fi

if [[ -f "$ANNOT_BED" ]]; then
  MATRIX_GZ="${BAMQC_DIR}/read_distribution_matrix.gz"
  PROFILE_PNG="${BAMQC_DIR}/read_distribution_profile.png"
  if [[ ! -f "$MATRIX_GZ" ]]; then
    computeMatrix scale-regions -S "${sorted_bams[@]}" -R "$ANNOT_BED" \
      --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
      -o "$MATRIX_GZ" -p "$THREADS"
  fi
  if [[ ! -f "$PROFILE_PNG" ]]; then
    plotProfile -m "$MATRIX_GZ" -out "$PROFILE_PNG"
  fi
fi

NPZ="${BAMQC_DIR}/multiBamSummary.npz"
if [[ -f "$NPZ" ]]; then
  labels=""
  for bam in "${sorted_bams[@]}"; do
    s=$(basename "$bam" _sorted.bam)
    labels+="$s "
  done
  PCA_PNG="${BAMQC_DIR}/deeptools_PCA.png"
  CORR_PNG="${BAMQC_DIR}/deeptools_correlation_heatmap.png"

  [[ ! -f "$PCA_PNG" ]] && plotPCA -in "$NPZ" --labels $labels -o "$PCA_PNG"
  [[ ! -f "$CORR_PNG" ]] && plotCorrelation -in "$NPZ" --corMethod spearman --skipZeros \
      --whatToPlot heatmap --colorMap RdYlBu_r --plotNumbers --labels $labels -o "$CORR_PNG"
fi

###############################################################################
# STEP 3: MACS2 summary + FRiP + FragLen
###############################################################################

echo "=== STEP 3: MACS2 summary + FRiP + Fragment Length ==="

# We rename the file slightly to avoid 'macs2' keyword collisions in auto-detection
MACS2_SUMMARY="${MACS2_DIR}/custom_macs2_stats_mqc.tsv"

# Header with Capitalized names
echo -e "Sample\tNum_Peaks\tPeak_Type\tFRiP\tFrag_Len" > "$MACS2_SUMMARY"

shopt -s nullglob
for bed in "${MACS2_DIR}"/*_filtered.bed; do
  name=$(basename "$bed" _filtered.bed)

  # 1. Peak Type
  peak_type="Unknown"
  [[ -f "${MACS2_DIR}/${name}_peaks.narrowPeak" ]] && peak_type="Narrow"
  [[ -f "${MACS2_DIR}/${name}_peaks.broadPeak" ]] && peak_type="Broad"

  # 2. Number of Peaks
  n_peaks=$(wc -l < "$bed")

  # 3. FRiP (Read from UPSTREAM generated table)
  frip="NA"
  if [[ -f "$FRIP_TSV" ]]; then
    # The upstream script writes: name \t reads_in_peaks \t total_mapped \t frip
    frip=$(awk -v s="$name" '$1==s {print $4}' "$FRIP_TSV")
    [[ -z "$frip" ]] && frip="NA"
  fi

  # 4. Fragment Length (from _peaks.xls comment "# d = ...")
  frag_len="NA"
  xls_file="${MACS2_DIR}/${name}_peaks.xls"
  if [[ -f "$xls_file" ]]; then
    d_val=$(grep -m 1 "# d = " "$xls_file" | sed 's/.*# d = //')
    [[ -n "$d_val" ]] && frag_len="$d_val"
  fi

  echo -e "${name}\t${n_peaks}\t${peak_type}\t${frip}\t${frag_len}" >> "$MACS2_SUMMARY"
done
shopt -u nullglob

###############################################################################
# STEP 4: HOMER annotatePeaks.pl
###############################################################################

echo "=== STEP 4: HOMER annotatePeaks + category summary ==="

mkdir -p "$HOMER_ANN_DIR"

# (Re-running generation loop just in case, logic unchanged)
shopt -s nullglob
for bed in "${MACS2_DIR}"/*_filtered.bed; do
  name=$(basename "$bed" _filtered.bed)
  ann_out="${HOMER_ANN_DIR}/${name}.annotatePeaks.txt"

  if [[ ! -f "$ann_out" ]]; then
    echo "[STEP4] Annotating $name..."
    annotatePeaks.pl "$bed" "$GENOME" > "$ann_out"
  fi
done
shopt -u nullglob

HOMER_SUMMARY="${HOMER_ANN_DIR}/homer_annotation_summary_mqc.tsv"
echo -e "Sample\tPromoter_TSS\tUTR5\tExon\tIntron\tIntergenic\tTTS" > "$HOMER_SUMMARY"

shopt -s nullglob
for ann in "${HOMER_ANN_DIR}"/*.annotatePeaks.txt; do
  sample=$(basename "$ann" .annotatePeaks.txt)
  
  awk -v s="$sample" '
    BEGIN { FS="\t"; OFS="\t"; prom=0; utr5=0; exon=0; intron=0; intergenic=0; tts=0; total=0; }
    NR==1 { next }
    {
      total++;
      ann=$8;
      if (ann ~ /Promoter|TSS/)       prom++;
      else if (ann ~ /5.?UTR/)        utr5++;
      else if (ann ~ /Exon/)          exon++;
      else if (ann ~ /Intron/)        intron++;
      else if (ann ~ /TTS/)           tts++;
      else                            intergenic++;
    }
    END {
      if(total==0) total=1;
      printf "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", s, 
      100*prom/total, 100*utr5/total, 100*exon/total, 
      100*intron/total, 100*intergenic/total, 100*tts/total;
    }
  ' "$ann" >> "$HOMER_SUMMARY"
done
shopt -u nullglob

###############################################################################
# STEP 5: Write multiqc_config.yaml
###############################################################################

echo "=== STEP 5: Writing multiqc_config.yaml ==="

cat > "${BASE}/multiqc_config.yaml" <<EOF
# 1. EXCLUDE MODULES
exclude_modules:
  - homer

# 2. REPORT SECTION ORDER (Vertical order of the report graphs)
#    Higher numbers = Higher up on the page
report_section_order:
  fastqc:
    order: 1000
  cutadapt:
    order: 900
  # REQUEST: Bowtie2 before Samtools
  bowtie2:
    order: 850
  samtools:
    order: 800
  # Standard MACS2 plots
  macs:
    order: 750
  # Our Custom Table (Detail view with FRiP)
  custom_macs2_stats:
    order: 700
  phantompeakqualtools:
    order: 600
  homer_annotation:
    order: 500

# 3. TABLE COLUMN PLACEMENT (Left-to-Right order)
#    Low numbers = Left side
table_columns_placement:
  fastqc:
    total_sequences: 10
    percent_gc: 20
    percent_duplicates: 30
  cutadapt:
    percent_trimmed: 40
  samtools_stats:
    reads_mapped_percent: 50
  # Standard MACS2 (Native module)
  macs2:
    peak_count: 60
    d: 70
  # SPP
  phantompeakqualtools:
    NSC: 80
    RSC: 90

# 4. VISIBILITY RULES (The Whitelist)
table_columns_visible:
  # --- FastQC ---
  fastqc-total_sequences: True
  fastqc-percent_gc: True
  fastqc-percent_duplicates: True
  
  # --- Cutadapt ---
  cutadapt-percent_trimmed: True
  
  # --- Samtools ---
  samtools_stats-reads_mapped_percent: True
  
  # --- MACS2 (Standard Module) ---
  # REQUEST: Only show Fragment Length (d) and Peak Count
  macs2-peak_count: True
  macs2-d: True
  
  # HIDE these MACS2 columns
  macs2-treatment_redundant_rate: False
  macs2-control_redundant_rate: False
  macs2-num_peaks: False   # Just in case of variable naming conflicts

  # --- SPP ---
  phantompeakqualtools-NSC: True
  phantompeakqualtools-RSC: True
  phantompeakqualtools-Estimated_Fragment_Length_bp: False # Hiding this as per request
  
  # --- Hide Clutter ---
  samtools_stats-error_rate: False
  samtools_stats-non_primary_alignments: False
  samtools_stats-reads_mapped: False
  samtools_stats-reads_properly_paired_percent: False
  samtools_stats-raw_total_sequences: False
  samtools_stats-insert_size_average: False
  samtools_flagstat-flagstat_total: False
  samtools_flagstat-mapped_passed: False
  samtools_flagstat-mapped_passed_pct: False
  bowtie_2_hisat2-overall_alignment_rate: False

# 5. COLUMN FORMATTING
general_stats_table_table:
  - fastqc:
      total_sequences:
        name: "M Seqs"
        scale: 1e-6
        format: "{:,.1f}"
      percent_gc:
        name: "%GC"
      percent_duplicates:
        name: "%Dups"
  
  - cutadapt:
      percent_trimmed:
        name: "%Trim"

  - samtools_stats:
      reads_mapped_percent:
        name: "%Mapped"
        format: "{:.1f}%"

  - macs2:
      peak_count:
        name: "MACS2 Peaks"
        format: "{:,.0f}"
      d:
        name: "Frag Len"
        suffix: " bp"
        format: "{:.0f}"

  - phantompeakqualtools:
      NSC:
        name: "NSC"
        format: "{:.2f}"
      RSC:
        name: "RSC"
        format: "{:.2f}"

# 6. CUSTOM DATA (Still keeping the Custom MACS2 table for the Detail section)
custom_data:
  custom_macs2_stats:
    file: "peaks/macs2/custom_macs2_stats_mqc.tsv"
    section_name: "MACS2 Peak Summary"
    description: "Detailed peak stats including FRiP (Fraction of Reads in Peaks)"
    plot_type: "table"
    pconfig:
      id: "custom_macs2_stats"
      title: "MACS2 Peak Metrics"
  
  homer_annotation:
    file: "peaks/homer/homer_annotation_summary_mqc.tsv"
    section_name: "HOMER Peak Annotation"
    description: "Peak annotation proportions"
    plot_type: "bargraph"
    pconfig:
      id: "homer_annotation"
      title: "Peak Annotation Categories"
      ylab: "Percent Peaks"
      stacking: "normal"

sample_names_replace_regex: true
sample_names_replace:
  "_R[12](_val_[12])?$": ""
  "_sorted$": ""
  "_bowtie2$": ""
  "_vs_.*$": ""

fn_ignore_files:
  - "*.sh"
  - "*multiqc_config.yaml"
EOF

###############################################################################
# STEP 6: Run MultiQC
###############################################################################

echo "=== STEP 6: Running MultiQC ==="

multiqc "$BASE" -o "$MULTIQC_OUT" -c "${BASE}/multiqc_config.yaml" -f

echo "=== CUT&RUN extra QC + MultiQC finished successfully ==="