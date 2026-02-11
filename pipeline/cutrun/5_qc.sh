#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 5_qc.sh — CUT&RUN QC + MultiQC
#
# Run after alignment (steps 1–4). Requires sorted BAMs in BAM_DIR.
#
# Steps:
#   1) Per-BAM: preseq (complexity extrapolation) + phantompeakqualtools (run_spp.R)
#   2) deepTools: plotFingerprint, multiBamSummary → plotPCA/plotCorrelation
#   3) MACS3 summary (optional): only if MACS3_DIR has peak files
#   4) Write multiqc_config.yaml (minimal; merges sample names)
#   5) Run multiqc
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "5_multiqc"

for c in samtools preseq run_spp.R plotFingerprint plotPCA plotCorrelation multiBamSummary multiqc; do
  check_cmd "$c"
done

THREADS="${QC_THREADS:-8}"
# BAMQC_DIR and MULTIQC_OUT are created in 0_config.sh; mkdir -p here is redundant but harmless

echo "=== CUT&RUN QC + MultiQC ==="

###############################################################################
# Inputs: sorted BAMs
###############################################################################

sorted_bams=( "${BAM_DIR}"/*_sorted.bam )

if (( ${#sorted_bams[@]} == 0 )); then
  echo "[WARN] No *_sorted.bam in $BAM_DIR. Exiting."
  exit 0
fi

###############################################################################
# STEP 1: Per-BAM QC — preseq + phantompeakqualtools
###############################################################################

echo "=== STEP 1: preseq + phantompeakqualtools ==="

for bam in "${sorted_bams[@]}"; do
  sample=$(basename "$bam" _sorted.bam)
  echo "[STEP1] $sample"

  preseq_file="${BAMQC_DIR}/${sample}.preseq.txt"
  if [[ ! -f "$preseq_file" ]]; then
    echo "  - preseq lc_extrap -> $preseq_file"
    preseq lc_extrap -B "$bam" -o "$preseq_file" -v
  fi

  spp_out="${BAMQC_DIR}/${sample}.spp.out"
  spp_pdf="${BAMQC_DIR}/${sample}.spp.pdf"
  if [[ ! -f "$spp_out" ]]; then
    echo "  - run_spp.R -> $spp_out"
    run_spp.R -c="$bam" -p="$THREADS" -savp="$spp_pdf" -out="$spp_out"
  fi
done

###############################################################################
# STEP 2: deepTools — plotFingerprint, multiBamSummary, PCA/Correlation
###############################################################################

echo "=== STEP 2: deepTools ==="

# plotFingerprint
FINGERPRINT_RAW="${BAMQC_DIR}/fingerprint_rawcounts.txt"
FINGERPRINT_PNG="${BAMQC_DIR}/fingerprint.png"
if [[ ! -f "$FINGERPRINT_RAW" ]]; then
  echo "[STEP2] plotFingerprint"
  plotFingerprint -b "${sorted_bams[@]}" -o "$FINGERPRINT_PNG" \
    --outRawCounts "$FINGERPRINT_RAW" --ignoreDuplicates -p "$THREADS"
fi

# multiBamSummary (for plotPCA / plotCorrelation)
NPZ="${BAMQC_DIR}/multiBamSummary.npz"
if [[ ! -f "$NPZ" ]]; then
  echo "[STEP2] multiBamSummary -> $NPZ"
  if [[ -n "${BAMQC_REGION:-}" ]]; then
    multiBamSummary bins -b "${sorted_bams[@]}" --region "$BAMQC_REGION" --outFileName "$NPZ" -p "$THREADS"
  else
    multiBamSummary bins -b "${sorted_bams[@]}" --outFileName "$NPZ" -p "$THREADS"
  fi
fi

# plotPCA / plotCorrelation (output data files so MultiQC can draw PCA and heatmap itself)
labels=()
for bam in "${sorted_bams[@]}"; do
  labels+=( "$(basename "$bam" _sorted.bam)" )
done
PCA_PNG="${BAMQC_DIR}/deeptools_PCA.png"
PCA_DATA="${BAMQC_DIR}/deeptools_PCA.tab"
CORR_PNG="${BAMQC_DIR}/deeptools_correlation_heatmap.png"
CORR_MATRIX="${BAMQC_DIR}/deeptools_correlation_matrix.tab"
if [[ ! -f "$PCA_PNG" || ! -f "$PCA_DATA" ]]; then
  echo "[STEP2] plotPCA"
  plotPCA -in "$NPZ" --labels "${labels[@]}" -o "$PCA_PNG" --outFileNameData "$PCA_DATA"
fi
if [[ ! -f "$CORR_PNG" || ! -f "$CORR_MATRIX" ]]; then
  echo "[STEP2] plotCorrelation"
  plotCorrelation -in "$NPZ" --corMethod spearman --skipZeros \
    --whatToPlot heatmap --plotNumbers --labels "${labels[@]}" -o "$CORR_PNG" \
    --outFileCorMatrix "$CORR_MATRIX"
fi

###############################################################################
# STEP 3: MACS3 summary (optional — only if MACS3_DIR has content)
###############################################################################

MACS_SUMMARY="${MACS3_DIR}/custom_macs_stats_mqc.tsv"
macs_has_peaks=0
for bed in "${MACS3_DIR}"/*_filtered.bed; do
  [[ -f "$bed" ]] && { macs_has_peaks=1; break; }
done

if [[ "$macs_has_peaks" -eq 1 ]]; then
  echo "=== STEP 3: MACS3 summary ==="
  echo -e "Sample\tNum_Peaks\tPeak_Type\tFRiP\tFrag_Len" > "$MACS_SUMMARY"
  for bed in "${MACS3_DIR}"/*_filtered.bed; do
    [[ -f "$bed" ]] || continue
    name=$(basename "$bed" _filtered.bed)

    peak_type="Unknown"
    [[ -f "${MACS3_DIR}/${name}_peaks.narrowPeak" ]] && peak_type="Narrow"
    [[ -f "${MACS3_DIR}/${name}_peaks.broadPeak"  ]] && peak_type="Broad"

    n_peaks=$(wc -l < "$bed" | tr -d ' ')

    frag_len="NA"
    xls_file="${MACS3_DIR}/${name}_peaks.xls"
    if [[ -f "$xls_file" ]]; then
      d_val=$(grep -m 1 "# d = " "$xls_file" 2>/dev/null | sed 's/.*# d = //')
      [[ -n "${d_val:-}" ]] && frag_len="$d_val"
    fi

    frip="NA"
    ip_bam="${BAM_DIR}/${name}_sorted.bam"
    if [[ -f "$ip_bam" ]]; then
      if [[ -f "${ip_bam}.bai" || -f "${ip_bam%.bam}.bai" ]]; then
        total_mapped=$(samtools view -c -F 260 "$ip_bam")
        reads_in_peaks=$(samtools view -c -F 260 -L "$bed" "$ip_bam")
        frip=$(awk -v inpeaks="$reads_in_peaks" -v total="$total_mapped" \
          'BEGIN { if (total > 0) printf "%.4f", inpeaks/total; else print "0.0000" }')
      fi
    fi

    echo -e "${name}\t${n_peaks}\t${peak_type}\t${frip}\t${frag_len}" >> "$MACS_SUMMARY"
  done
else
  echo "=== STEP 3: MACS3 dir empty or no *_filtered.bed; skipping MACS3 summary ==="
fi

###############################################################################
# STEP 4: multiqc_config.yaml (sample names, General Statistics columns)
###############################################################################

echo "=== STEP 4: Writing multiqc_config.yaml ==="

# Build config: base + optional MACS custom_data
MACS_CUSTOM=""
if [[ -f "$MACS_SUMMARY" && -s "$MACS_SUMMARY" ]]; then
  MACS_CUSTOM="
# MACS3 custom table (optional)
custom_data:
  custom_macs_stats:
    file: \"peaks/macs3/custom_macs_stats_mqc.tsv\"
    section_name: \"MACS Peak Summary\"
    plot_type: \"table\"
    pconfig:
      id: \"custom_macs_stats\"
      title: \"MACS Peak Metrics\"
"
fi

cat > "${BASE}/multiqc_config.yaml" <<EOF
# Merge sample names so data, data_R1, data_R1_val_1, data_R2, data_R2_val_2, data_sorted, etc. -> one sample
sample_names_replace_regex: true
sample_names_replace:
  '_R[12](_val_[12])?\\.fq\\.gz$': ""
  '_R[12](_val_[12])?$': ""
  '_val_[12]$': ""
  '_sorted$': ""
  '_bowtie2$': ""
  '\\.bam$': ""
  '_vs_.*\$': ""

fn_ignore_files:
  - "*.sh"
  - "*multiqc_config.yaml"
  - "*.bam"
  - "*.bai"

# Module/section order in report (top to bottom). Lower order value = appears earlier.
# Find section IDs by clicking navigation links in the report (URL shows #section_id).
report_section_order:
  general_stats:
    order: 1000
  fastqc:
    order: 900
  cutadapt:
    order: 850
  bowtie2:
    order: 800
  samtools:
    order: 750
  deeptools:
    order: 700
  homer:
    order: 600
  phantompeakqualtools:
    order: 500
  preseq:
    order: 400
  custom_content:
    order: 300

# General Statistics: column order (left to right). Lower placement value = left.
# Format: module namespace (as shown in Configure Columns) -> column ID -> placement value.
table_columns_placement:
  fastqc:
    total_sequences: 100
    percent_duplicates: 200
    percent_gc: 300
  cutadapt:
    percent_trimmed: 400
  bowtie_2_hisat2:
    overall_alignment_rate: 500
  homer:
    averageTagLength: 600
    peakSizeEstimate: 700
    fragmentLengthEstimate: 800
  homer_findpeaks:
    total_peaks: 900
    approximate_ip_efficiency: 1000
  phantompeakqualtools:
    NSC: 1100
    RSC: 1200

# Show only these columns; hide all others (including entire samtools_stats and unwanted homer/SPP columns).
table_columns_visible:
  fastqc:
    total_sequences: true
    percent_duplicates: true
    percent_gc: true
  cutadapt:
    percent_trimmed: true
  bowtie_2_hisat2:
    overall_alignment_rate: true
  homer:
    averageTagLength: true
    peakSizeEstimate: true
    fragmentLengthEstimate: true
    tagsPerBP: false
    averageFragmentGCcontent: false
  homer_findpeaks:
    total_peaks: true
    approximate_ip_efficiency: true
    expected_tags_per_peak: false
  phantompeakqualtools:
    NSC: true
    RSC: true
    Estimated_Fragment_Length_bp: false
  samtools_stats: false
  'Samtools: stats': false
  # Explicitly hide individual samtools_stats columns (fallback if namespace hiding doesn't work)
  samtools_stats-reads_mapped_percent: false
  samtools_stats-reads_properly_paired_percent: false
  samtools_stats-reads_mapped: false
  samtools_stats-error_rate: false
  samtools_stats-non_primary_alignments: false
  samtools_stats-reads_MQ0_percent: false
  samtools_stats-insert_size_average: false
  samtools_stats-raw_total_sequences: false
  # Hide homer_findpeaks-expected_tags_per_peak (try flat format if nested doesn't work)
  'homer_findpeaks-expected_tags_per_peak': false
${MACS_CUSTOM}
EOF

###############################################################################
# STEP 5: Run MultiQC (BASE + BAMQC_DIR; deepTools PCA/correlation from .tab data files)
###############################################################################

echo "=== STEP 5: Running MultiQC ==="
multiqc "$BASE" "$BAMQC_DIR" -o "$MULTIQC_OUT" -c "${BASE}/multiqc_config.yaml" -f

echo "=== 5_qc finished successfully ==="
