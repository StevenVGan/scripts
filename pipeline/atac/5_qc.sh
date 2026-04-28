#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 5_qc.sh — ATAC QC + MultiQC
#
# Differences vs cutrun/5_qc.sh:
#   - Drops run_spp.R (designed for ChIP, not informative for ATAC)
#   - Adds bamPEFragmentSize (NFR/mono/di/tri-nucleosomal periodicity check)
#   - Adds TSS enrichment (computeMatrix reference-point on bigWigs +/- TSS_FLANK)
#   - Auto-generates the TSS BED from $HOME/work/ref/annot/${GENOME}/${GENOME}_annot_genome.bed
#     if $TSS_BED is missing
#   - preseq runs on ${sample}_marked.bam (pre-dedup) so complexity extrapolation
#     reflects the full library; falls back to ${sample}_sorted.bam if marked is gone
#   - MACS3 summary now expects narrowPeak only (no broad branch)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "5_qc"

for c in samtools preseq plotFingerprint plotPCA plotCorrelation multiBamSummary \
         bamPEFragmentSize computeMatrix plotProfile multiqc; do
  check_cmd "$c"
done

THREADS="${QC_THREADS:-8}"

echo "=== ATAC QC + MultiQC ==="

###############################################################################
# Inputs
###############################################################################

sorted_bams=( "${BAM_DIR}"/*_sorted.bam )

if (( ${#sorted_bams[@]} == 0 )); then
  echo "[WARN] No *_sorted.bam in $BAM_DIR. Exiting."
  exit 0
fi

###############################################################################
# STEP 1: preseq (per-BAM library complexity)
#         Prefer _marked.bam (pre-dedup) for accurate extrapolation
###############################################################################

echo "=== STEP 1: preseq ==="

for bam in "${sorted_bams[@]}"; do
  sample=$(basename "$bam" _sorted.bam)

  preseq_in="${BAM_DIR}/${sample}_marked.bam"
  if [[ ! -f "$preseq_in" ]]; then
    echo "[STEP1] $sample: marked.bam missing, falling back to sorted.bam (will underestimate dup rate)"
    preseq_in="$bam"
  fi

  preseq_file="${BAMQC_DIR}/${sample}.preseq.txt"
  if [[ ! -f "$preseq_file" ]]; then
    echo "[STEP1] $sample preseq lc_extrap -> $preseq_file"
    preseq lc_extrap -B "$preseq_in" -o "$preseq_file" -v
  fi
done

###############################################################################
# STEP 2: deepTools — fingerprint, PCA/correlation
###############################################################################

echo "=== STEP 2: deepTools fingerprint + PCA/correlation ==="

FINGERPRINT_RAW="${BAMQC_DIR}/fingerprint_rawcounts.txt"
FINGERPRINT_PNG="${BAMQC_DIR}/fingerprint.png"
if [[ ! -f "$FINGERPRINT_RAW" ]]; then
  plotFingerprint -b "${sorted_bams[@]}" -o "$FINGERPRINT_PNG" \
    --outRawCounts "$FINGERPRINT_RAW" --ignoreDuplicates -p "$THREADS"
fi

NPZ="${BAMQC_DIR}/multiBamSummary.npz"
if [[ ! -f "$NPZ" ]]; then
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
  plotPCA -in "$NPZ" --labels "${labels[@]}" -o "$PCA_PNG" --outFileNameData "$PCA_DATA"
fi
if [[ ! -f "$CORR_PNG" || ! -f "$CORR_MATRIX" ]]; then
  plotCorrelation -in "$NPZ" --corMethod spearman --skipZeros \
    --whatToPlot heatmap --plotNumbers --labels "${labels[@]}" -o "$CORR_PNG" \
    --outFileCorMatrix "$CORR_MATRIX"
fi

###############################################################################
# STEP 3: ATAC fragment-size distribution (NFR / mononuc / dinuc / trinuc)
###############################################################################

if [[ "${SE:-0}" -eq 0 ]]; then
  echo "=== STEP 3: bamPEFragmentSize ==="
  FRAG_PNG="${BAMQC_DIR}/fragment_size.png"
  FRAG_TSV="${BAMQC_DIR}/fragment_size.tsv"
  FRAG_RAW="${BAMQC_DIR}/fragment_size_raw.tsv"
  if [[ ! -f "$FRAG_PNG" ]]; then
    bamPEFragmentSize -b "${sorted_bams[@]}" --samplesLabel "${labels[@]}" \
      -hist "$FRAG_PNG" --table "$FRAG_TSV" --outRawFragmentLengths "$FRAG_RAW" \
      -p "$THREADS"
  fi
else
  echo "=== STEP 3: SE mode — skipping bamPEFragmentSize ==="
fi

###############################################################################
# STEP 4: TSS enrichment (auto-generate TSS BED if missing)
###############################################################################

echo "=== STEP 4: TSS enrichment ==="

if [[ ! -f "$TSS_BED" ]]; then
  src_bed="${HOME}/work/ref/annot/${GENOME}/${GENOME}_annot_genome.bed"
  if [[ -f "$src_bed" ]]; then
    echo "[STEP4] Generating TSS BED from $src_bed -> $TSS_BED"
    mkdir -p "$(dirname "$TSS_BED")"
    # Per-transcript TSS: + strand uses BED start, - strand uses BED end-1
    awk -v OFS='\t' '
      $6 == "+" { print $1, $2,   $2+1, $4, $5, $6 }
      $6 == "-" { print $1, $3-1, $3,   $4, $5, $6 }
    ' "$src_bed" | sort -k1,1 -k2,2n -u > "$TSS_BED"
  else
    echo "[WARN] No source annotation at $src_bed; skipping TSS enrichment."
    TSS_BED=""
  fi
fi

bigwigs=( "${TRACK_DIR}"/*.bw )
if [[ -n "$TSS_BED" && -f "$TSS_BED" && ${#bigwigs[@]} -gt 0 ]]; then
  TSS_MATRIX="${BAMQC_DIR}/tss_matrix.gz"
  TSS_PROFILE="${BAMQC_DIR}/tss_profile.png"
  TSS_TABLE="${BAMQC_DIR}/tss_profile_data.tab"
  bw_labels=()
  for bw in "${bigwigs[@]}"; do bw_labels+=( "$(basename "$bw" .bw)" ); done

  if [[ ! -f "$TSS_MATRIX" ]]; then
    computeMatrix reference-point -S "${bigwigs[@]}" -R "$TSS_BED" \
      --referencePoint center -a "$TSS_FLANK" -b "$TSS_FLANK" \
      --skipZeros -p "$THREADS" -o "$TSS_MATRIX"
  fi
  if [[ ! -f "$TSS_PROFILE" ]]; then
    plotProfile -m "$TSS_MATRIX" -o "$TSS_PROFILE" \
      --samplesLabel "${bw_labels[@]}" --outFileNameData "$TSS_TABLE" \
      --regionsLabel "TSS" --plotTitle "ATAC signal at TSS"
  fi
fi

###############################################################################
# STEP 5: MACS3 summary (narrow only; FRiP from filtered _sorted.bam)
###############################################################################

MACS_SUMMARY="${MACS3_DIR}/custom_macs_stats_mqc.tsv"
macs_has_peaks=0
for bed in "${MACS3_DIR}"/*_filtered.bed; do
  [[ -f "$bed" ]] && { macs_has_peaks=1; break; }
done

if [[ "$macs_has_peaks" -eq 1 ]]; then
  echo "=== STEP 5: MACS3 summary ==="
  echo -e "Sample\tNum_Peaks\tFRiP" > "$MACS_SUMMARY"
  for bed in "${MACS3_DIR}"/*_filtered.bed; do
    [[ -f "$bed" ]] || continue
    name=$(basename "$bed" _filtered.bed)
    n_peaks=$(wc -l < "$bed" | tr -d ' ')

    frip="NA"
    ip_bam="${BAM_DIR}/${name}_sorted.bam"
    if [[ -f "$ip_bam" && ( -f "${ip_bam}.bai" || -f "${ip_bam%.bam}.bai" ) ]]; then
      total_mapped=$(samtools view -c -F 260 "$ip_bam")
      reads_in_peaks=$(samtools view -c -F 260 -L "$bed" "$ip_bam")
      frip=$(awk -v inpeaks="$reads_in_peaks" -v total="$total_mapped" \
        'BEGIN { if (total > 0) printf "%.4f", inpeaks/total; else print "0.0000" }')
    fi

    echo -e "${name}\t${n_peaks}\t${frip}" >> "$MACS_SUMMARY"
  done
else
  echo "=== STEP 5: MACS3 dir empty; skipping summary ==="
fi

###############################################################################
# STEP 6: multiqc_config.yaml
###############################################################################

echo "=== STEP 6: Writing multiqc_config.yaml ==="

MACS_CUSTOM=""
if [[ -f "$MACS_SUMMARY" && -s "$MACS_SUMMARY" ]]; then
  MACS_CUSTOM="
custom_data:
  custom_macs_stats:
    file: \"peaks/macs3/custom_macs_stats_mqc.tsv\"
    section_name: \"MACS Peak Summary (ATAC)\"
    plot_type: \"table\"
    pconfig:
      id: \"custom_macs_stats\"
      title: \"MACS Peak Metrics\"
"
fi

cat > "${BASE}/multiqc_config.yaml" <<EOF
sample_names_replace_regex: true
sample_names_replace:
  '_R[12](_val_[12])?\\.fq\\.gz$': ""
  '_R[12](_val_[12])?$': ""
  '_val_[12]$': ""
  '_sorted\\.shifted$': ""
  '_sorted$': ""
  '_marked$': ""
  '_bowtie2$': ""
  '\\.bam$': ""

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
    order: 400
  custom_content:
    order: 300

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
    expected_tags_per_peak: false
  samtools_stats: false
  'Samtools: stats': false
${MACS_CUSTOM}
EOF

###############################################################################
# STEP 7: Run MultiQC
###############################################################################

echo "=== STEP 7: Running MultiQC ==="
multiqc "$BASE" "$BAMQC_DIR" -o "$MULTIQC_OUT" -c "${BASE}/multiqc_config.yaml" -f

echo "=== 5_qc finished successfully ==="
