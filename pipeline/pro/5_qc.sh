#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 5_qc.sh — PRO-seq QC + MultiQC
#
# Run after alignment (steps 1–4). Requires sorted BAMs in BAM_DIR.
#
# Steps:
#   1) Per-BAM: preseq (library-complexity extrapolation). phantompeakqualtools
#      is intentionally omitted — NSC/RSC use a ChIP-style fragment-bimodal
#      cross-correlation signature that does not fit TSS/nascent-RNA data.
#   2) deepTools: plotFingerprint (--ignoreDuplicates PE only), multiBamSummary
#      → plotPCA / plotCorrelation.
#   3) MACS3 summary (optional). Pausing-index summary (optional; produced by
#      step 4.3 as pausing_index_summary_mqc.tsv).
#   4) multiqc_config.yaml (sample-name merges, column order/visibility).
#   5) Run multiqc.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "5_multiqc"

for c in samtools preseq plotFingerprint plotPCA plotCorrelation multiBamSummary multiqc; do
  check_cmd "$c"
done

THREADS="${QC_THREADS:-8}"

echo "=== PRO-seq QC + MultiQC ==="

###############################################################################
# Inputs: sorted BAMs
###############################################################################

sorted_bams=( "${BAM_DIR}"/*_sorted.bam )

if (( ${#sorted_bams[@]} == 0 )); then
  echo "[WARN] No *_sorted.bam in $BAM_DIR. Exiting."
  exit 0
fi

###############################################################################
# STEP 1: preseq
###############################################################################

echo "=== STEP 1: preseq (library complexity) ==="

for bam in "${sorted_bams[@]}"; do
  sample=$(basename "$bam" _sorted.bam)
  echo "[STEP1] $sample"

  preseq_file="${BAMQC_DIR}/${sample}.preseq.txt"
  if [[ ! -f "$preseq_file" ]]; then
    echo "  - preseq lc_extrap -> $preseq_file"
    preseq lc_extrap -B "$bam" -o "$preseq_file" -v
  fi
done

###############################################################################
# STEP 2: deepTools — plotFingerprint, multiBamSummary, PCA/Correlation
###############################################################################

echo "=== STEP 2: deepTools ==="

fp_extra=()
[[ "${SE:-0}" -ne 1 ]] && fp_extra+=( --ignoreDuplicates )

FINGERPRINT_RAW="${BAMQC_DIR}/fingerprint_rawcounts.txt"
FINGERPRINT_PNG="${BAMQC_DIR}/fingerprint.png"
if [[ ! -f "$FINGERPRINT_RAW" ]]; then
  echo "[STEP2] plotFingerprint (SE=${SE:-0}${fp_extra[*]:+, }${fp_extra[*]:-})"
  plotFingerprint -b "${sorted_bams[@]}" -o "$FINGERPRINT_PNG" \
    --outRawCounts "$FINGERPRINT_RAW" ${fp_extra[@]+"${fp_extra[@]}"} -p "$THREADS"
fi

NPZ="${BAMQC_DIR}/multiBamSummary.npz"
if [[ ! -f "$NPZ" ]]; then
  echo "[STEP2] multiBamSummary -> $NPZ"
  if [[ -n "${BAMQC_REGION:-}" ]]; then
    multiBamSummary bins -b "${sorted_bams[@]}" --region "$BAMQC_REGION" --outFileName "$NPZ" -p "$THREADS"
  else
    multiBamSummary bins -b "${sorted_bams[@]}" --outFileName "$NPZ" -p "$THREADS"
  fi
fi

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
# STEP 3b: Pausing-index summary (optional — produced by step 4.3)
###############################################################################

PI_QUART_SUMMARY="${PAUSING_DIR}/pausing_index_quartiles_mqc.tsv"
PI_NGENES_SUMMARY="${PAUSING_DIR}/pausing_index_ngenes_mqc.tsv"
if [[ -s "$PI_QUART_SUMMARY" || -s "$PI_NGENES_SUMMARY" ]]; then
  echo "=== STEP 3b: Pausing-index summary present (quartiles + N_genes split) ==="
else
  echo "=== STEP 3b: No pausing_index_*_mqc.tsv (step 4.3 not run); skipping ==="
fi

###############################################################################
# STEP 3c: HOMER transcript summary (step 4.2 findPeaks -style groseq output)
#
# Built-in MultiQC homer_findpeaks module can miss -style groseq outputs, so we
# emit an explicit per-sample summary here. Columns:
#   Sample, N_transcripts (blacklist-filtered), total_peaks (raw findPeaks
#   header), approx_IP_efficiency (from findPeaks header).
###############################################################################

HOMER_SUMMARY="${HOMER_PEAK_DIR}/custom_homer_transcripts_mqc.tsv"
homer_has_output=0
for ann in "${HOMER_PEAK_DIR}"/*.annotatePeaks.txt; do
  [[ -f "$ann" ]] && { homer_has_output=1; break; }
done

if [[ "$homer_has_output" -eq 1 ]]; then
  echo "=== STEP 3c: HOMER transcript summary ==="
  echo -e "Sample\tN_transcripts\tTotal_peaks\tApprox_IP_efficiency" > "$HOMER_SUMMARY"
  for ann in "${HOMER_PEAK_DIR}"/*.annotatePeaks.txt; do
    [[ -f "$ann" ]] || continue
    name=$(basename "$ann" .annotatePeaks.txt)

    # Validated transcript count = non-header lines in blacklist-filtered annotatePeaks
    n_tx=$(tail -n +2 "$ann" | wc -l | tr -d ' ')

    # findPeaks headers preserved from 4.2 (either *_transcripts_stats.txt or others)
    stats_file=""
    for cand in "${HOMER_PEAK_DIR}/${name}_transcripts_stats.txt" \
                "${HOMER_PEAK_DIR}/${name}_peaks_stats.txt" \
                "${HOMER_PEAK_DIR}/${name}_regions_stats.txt"; do
      [[ -f "$cand" ]] && { stats_file="$cand"; break; }
    done

    total_peaks="NA"; ip_eff="NA"
    if [[ -n "$stats_file" ]]; then
      # findPeaks writes "# total peaks = N" (factor/histone) or "# total transcripts = N" (groseq).
      # `|| true` keeps set -euo pipefail from killing the script when the pattern is absent.
      tp=$(grep -m1 -E "^# total (peaks|transcripts) = " "$stats_file" 2>/dev/null | sed 's/.*= //' || true)
      [[ -n "$tp" ]] && total_peaks="$tp"
      # Approximate IP efficiency is only emitted when -i control is supplied; absent for groseq.
      eff=$(grep -m1 -iE "^# approximate IP efficiency = " "$stats_file" 2>/dev/null | sed 's/.*= //; s/%$//' || true)
      [[ -n "$eff" ]] && ip_eff="$eff"
    fi

    echo -e "${name}\t${n_tx}\t${total_peaks}\t${ip_eff}" >> "$HOMER_SUMMARY"
  done
else
  echo "=== STEP 3c: No HOMER annotatePeaks output; skipping transcript summary ==="
fi

###############################################################################
# STEP 4: multiqc_config.yaml
###############################################################################

echo "=== STEP 4: Writing multiqc_config.yaml ==="

CUSTOM_DATA_BLOCK=""
CUSTOM_HEADER="custom_data:"
if [[ -f "$MACS_SUMMARY" && -s "$MACS_SUMMARY" ]]; then
  CUSTOM_DATA_BLOCK+="
  custom_macs_stats:
    file: \"peaks/macs3/custom_macs_stats_mqc.tsv\"
    section_name: \"MACS Peak Summary\"
    plot_type: \"table\"
    pconfig:
      id: \"custom_macs_stats\"
      title: \"MACS Peak Metrics\"
"
fi
if [[ -f "$PI_QUART_SUMMARY" && -s "$PI_QUART_SUMMARY" ]]; then
  # PI quartile table: Median, Q1, Q3 share an axis (shared_key) so per-sample
  # bars are directly comparable across the three columns.
  CUSTOM_DATA_BLOCK+="
  pausing_index_quartiles:
    file: \"pausing/pausing_index_quartiles_mqc.tsv\"
    section_name: \"PRO-seq Pausing Index — quartiles\"
    description: \"Per-sample median, Q1 and Q3 of the per-gene pausing index (promoter density / gene-body density). Restricted to genes with ≥1 read in BOTH promoter and gene body. High median PI → pervasive Pol II promoter-proximal pausing.\"
    plot_type: \"table\"
    pconfig:
      id: \"pausing_index_quartiles\"
      title: \"Pausing Index — quartile summary\"
    headers:
      Median_PI:
        title: \"Median PI\"
        description: \"Median pausing index across genes\"
        format: \"{:,.2f}\"
        shared_key: \"pausing_index\"
      Q1_PI:
        title: \"Q1 PI\"
        format: \"{:,.2f}\"
        shared_key: \"pausing_index\"
      Q3_PI:
        title: \"Q3 PI\"
        format: \"{:,.2f}\"
        shared_key: \"pausing_index\"
"
fi
if [[ -f "$PI_NGENES_SUMMARY" && -s "$PI_NGENES_SUMMARY" ]]; then
  # Gene-count table is split off so its 10^5 magnitude doesn't crush the
  # PI quartile bars in a shared MultiQC violin/beeswarm view.
  CUSTOM_DATA_BLOCK+="
  pausing_index_ngenes:
    file: \"pausing/pausing_index_ngenes_mqc.tsv\"
    section_name: \"PRO-seq Pausing Index — gene count\"
    description: \"Number of genes contributing to the pausing-index summary above (≥1 read in BOTH promoter and gene body). Lower N_genes ≈ shallower depth or more aggressive zero-count filtering.\"
    plot_type: \"table\"
    pconfig:
      id: \"pausing_index_ngenes\"
      title: \"Genes evaluated for pausing index\"
    headers:
      N_genes:
        title: \"N genes\"
        description: \"Genes with ≥1 read in promoter AND gene body\"
        format: \"{:,.0f}\"
"
fi
if [[ -f "$HOMER_SUMMARY" && -s "$HOMER_SUMMARY" ]]; then
  CUSTOM_DATA_BLOCK+="
  custom_homer_transcripts:
    file: \"peaks/homer/custom_homer_transcripts_mqc.tsv\"
    section_name: \"HOMER findPeaks (step 4.2) — nascent transcripts\"
    description: \"Per-sample nascent-transcript counts from findPeaks -style ${HOMER_STYLE}; total_peaks and approx IP efficiency parsed from the findPeaks header block.\"
    plot_type: \"table\"
    pconfig:
      id: \"custom_homer_transcripts\"
      title: \"HOMER Transcript Metrics (findPeaks -style ${HOMER_STYLE})\"
    headers:
      N_transcripts:
        title: \"N transcripts\"
        description: \"Blacklist-filtered transcript count from annotatePeaks.txt\"
        format: \"{:,.0f}\"
      Total_peaks:
        title: \"Raw peaks\"
        description: \"# total peaks from findPeaks header (before annotation + blacklist filter)\"
        format: \"{:,.0f}\"
      Approx_IP_efficiency:
        title: \"IP efficiency %\"
        description: \"Approximate IP efficiency from findPeaks header\"
        format: \"{:,.2f}\"
"
fi
CUSTOM_BLOCK=""
if [[ -n "$CUSTOM_DATA_BLOCK" ]]; then
  CUSTOM_BLOCK="${CUSTOM_HEADER}${CUSTOM_DATA_BLOCK}"
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
  preseq:
    order: 500
  custom_content:
    order: 400

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
  samtools_stats: false
  'Samtools: stats': false
  samtools_stats-reads_mapped_percent: false
  samtools_stats-reads_properly_paired_percent: false
  samtools_stats-reads_mapped: false
  samtools_stats-error_rate: false
  samtools_stats-non_primary_alignments: false
  samtools_stats-reads_MQ0_percent: false
  samtools_stats-insert_size_average: false
  samtools_stats-raw_total_sequences: false
  'homer_findpeaks-expected_tags_per_peak': false
${CUSTOM_BLOCK}
EOF

###############################################################################
# STEP 5: Run MultiQC
###############################################################################

echo "=== STEP 5: Running MultiQC ==="
multiqc "$BASE" "$BAMQC_DIR" -o "$MULTIQC_OUT" -c "${BASE}/multiqc_config.yaml" -f

# Reference manifest — record genome/blacklist/bio_env_lock for reproducibility.
# See scripts/CONVENTIONS.md §4 + §9.
emit_references_tsv

echo "=== 5_qc (PRO-seq) finished successfully ==="
