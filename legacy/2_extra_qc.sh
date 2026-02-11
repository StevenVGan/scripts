#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# CUT&RUN EXTRA QC + MultiQC
#
# Run this AFTER your main pipeline has produced:
#   - align/*_sorted.bam (+ .bai)
#   - peaks/*_peaks.narrowPeak
#
# Generates:
#   - samtools flagstat/stats         → MultiQC alignment stats
#   - preseq lc_extrap                → library complexity curves
#   - deepTools plotFingerprint       → fingerprint plots
#   - phantompeakqualtools run_spp.R  → NSC/RSC/frag length
#   - HOMER annotatePeaks.pl          → peak annotation tables
#   - multiqc_config.yaml             → controls module order
#   - MultiQC report over the project
#
# Assumes conda env "bio" with:
#   samtools, preseq, deeptools, phantompeakqualtools, homer, multiqc
###############################################################################

# >>> EDIT THIS IF PROJECT PATH CHANGES <<<
BASE="${HOME}/work/seq/CUTRUN/250822_CnR_ERa_MCF7_delEZH2_Steven"

ALIGN_DIR="${BASE}/align"
PEAK_DIR="${BASE}/peaks"

# QC directories
QC_DIR="${ALIGN_DIR}/qc"
SAMTOOLS_QC_DIR="${QC_DIR}/samtools"
PRESEQ_QC_DIR="${QC_DIR}/preseq"
DEEPP_QC_DIR="${QC_DIR}/deeptools"
SPP_QC_DIR="${QC_DIR}/spp"
HOMER_QC_DIR="${PEAK_DIR}/homer"

# HOMER genome
HOMER_GENOME="hg38"

# Toggles (set 0 → skip step)
RUN_SAMTOOLS_QC=1
RUN_PRESEQ=1
RUN_FINGERPRINT=1
RUN_SPP_QC=1
RUN_HOMER_ANNOT=1
RUN_MULTIQC=1

# Threads / params
SPP_THREADS=8
FINGERPRINT_READS=500000

mkdir -p "$SAMTOOLS_QC_DIR" "$PRESEQ_QC_DIR" "$DEEPP_QC_DIR" "$SPP_QC_DIR" "$HOMER_QC_DIR"

# Helper to check commands
check_cmd() {
  local c="$1"
  if ! command -v "$c" >/dev/null 2>&1; then
    echo "ERROR: required command '$c' not found in PATH." >&2
    exit 1
  fi
}

echo "=== EXTRA QC starting (BASE=$BASE) ==="

# Check tools
check_cmd samtools
check_cmd preseq
check_cmd plotFingerprint
check_cmd run_spp.R
check_cmd annotatePeaks.pl
check_cmd multiqc

###############################################################################
# 1) Per-BAM QC
###############################################################################
echo "[INFO] Collecting BAMs from: $ALIGN_DIR"
sorted_bams=( "$ALIGN_DIR"/*_sorted.bam )

if (( ${#sorted_bams[@]} == 0 )); then
  echo "[WARN] No *_sorted.bam files found in $ALIGN_DIR; skipping BAM-based QC."
else
  for bam in "${sorted_bams[@]}"; do
    [[ -e "$bam" ]] || continue
    sample=$(basename "$bam" _sorted.bam)
    echo ">>> Sample: $sample"

    # 1.1 samtools flagstat / stats
    if [[ "$RUN_SAMTOOLS_QC" -eq 1 ]]; then
      flagstat_out="${SAMTOOLS_QC_DIR}/${sample}.flagstat.txt"
      stats_out="${SAMTOOLS_QC_DIR}/${sample}.stats.txt"

      if [[ ! -f "$flagstat_out" ]]; then
        echo "  [samtools flagstat]"
        samtools flagstat "$bam" > "$flagstat_out"
      else
        echo "  [samtools flagstat] exists, skipping"
      fi

      if [[ ! -f "$stats_out" ]]; then
        echo "  [samtools stats]"
        samtools stats "$bam" > "$stats_out"
      else
        echo "  [samtools stats] exists, skipping"
      fi
    fi

    # 1.2 preseq
    if [[ "$RUN_PRESEQ" -eq 1 ]]; then
      preseq_out="${PRESEQ_QC_DIR}/${sample}.preseq.dat"
      if [[ ! -f "$preseq_out" ]]; then
        echo "  [preseq lc_extrap]"
        preseq lc_extrap -B "$bam" -o "$preseq_out"
      else
        echo "  [preseq] exists, skipping"
      fi
    fi

    # 1.3 deepTools plotFingerprint
    if [[ "$RUN_FINGERPRINT" -eq 1 ]]; then
      fp_png="${DEEPP_QC_DIR}/${sample}_fingerprint.png"
      fp_raw="${DEEPP_QC_DIR}/${sample}_fingerprint.raw.txt"
      fp_metrics="${DEEPP_QC_DIR}/${sample}_fingerprint.metrics.txt"
      if [[ ! -f "$fp_png" ]]; then
        echo "  [deepTools plotFingerprint]"
        plotFingerprint \
          -b "$bam" \
          -o "$fp_png" \
          --outRawCounts "$fp_raw" \
          --outQualityMetrics "$fp_metrics" \
          --numberOfSamples "$FINGERPRINT_READS" \
          --skipZeros \
          || echo "  [WARN] plotFingerprint failed for $sample"
      else
        echo "  [plotFingerprint] exists, skipping"
      fi
    fi

    # 1.4 phantompeakqualtools (run_spp.R)
    if [[ "$RUN_SPP_QC" -eq 1 ]]; then
      spp_out="${SPP_QC_DIR}/${sample}.spp.out"
      spp_pdf="${SPP_QC_DIR}/${sample}.spp.pdf"
      if [[ ! -f "$spp_out" ]]; then
        echo "  [phantompeakqualtools run_spp.R]"
        run_spp.R \
          -c="$bam" \
          -p="$SPP_THREADS" \
          -savp \
          -out="$spp_out" \
          -savp="$spp_pdf" \
          || echo "  [WARN] run_spp.R failed for $sample"
      else
        echo "  [run_spp.R] exists, skipping"
      fi
    fi

  done
fi

###############################################################################
# 2) HOMER annotatePeaks on MACS2 peaks
###############################################################################
if [[ "$RUN_HOMER_ANNOT" -eq 1 ]]; then
  echo "=== HOMER annotatePeaks on MACS2 peaks ==="
  peak_files=( "$PEAK_DIR"/*_peaks.narrowPeak )

  if (( ${#peak_files[@]} == 0 )); then
    echo "[WARN] No *_peaks.narrowPeak files in $PEAK_DIR; skipping HOMER annotation."
  else
    for pf in "${peak_files[@]}"; do
      [[ -e "$pf" ]] || continue
      base=$(basename "$pf" _peaks.narrowPeak)
      out_annot="${HOMER_QC_DIR}/${base}.homer.annot.txt"

      if [[ ! -f "$out_annot" ]]; then
        echo "  [annotatePeaks.pl] $base"
        annotatePeaks.pl "$pf" "$HOMER_GENOME" > "$out_annot"
      else
        echo "  [annotatePeaks.pl] $base exists, skipping"
      fi
    done
  fi
fi

# --- Summarise HOMER annotatePeaks.pl output ---

HOMER_SUMMARY="${PEAK_DIR}/homer_peak_annotation.tsv"

# Header: sample + broad genomic categories (percent of peaks)
echo -e "sample\tPromoter\tExon\tIntron\tIntergenic\tDownstream\tUTR5\tUTR3\tTTS\tOther" > "$HOMER_SUMMARY"

for annot in "${PEAK_DIR}"/*_annot.txt; do
    [[ ! -f "$annot" ]] && continue
    sample=$(basename "$annot" _annot.txt)

    row=$(
        awk '
            BEGIN {
                FS = OFS = "\t"
            }
            NR==1 {
                # Find "Annotation" column
                for (i=1; i<=NF; i++) {
                    if ($i == "Annotation") {
                        ann_col = i
                        break
                    }
                }
                next
            }
            NR>1 {
                ann = $ann_col
                cat = "Other"
                if (ann ~ /Promoter/)         cat = "Promoter"
                else if (ann ~ /Exon/)        cat = "Exon"
                else if (ann ~ /Intron/)      cat = "Intron"
                else if (ann ~ /Intergenic/)  cat = "Intergenic"
                else if (ann ~ /Downstream/)  cat = "Downstream"
                else if (ann ~ /5.?UTR/)      cat = "UTR5"
                else if (ann ~ /3.?UTR/)      cat = "UTR3"
                else if (ann ~ /TTS/)         cat = "TTS"

                counts[cat]++
                total++
            }
            END {
                if (total == 0) total = 1
                # helper to get percentage or 0
                function pct(name) {
                    return (name in counts ? 100.0 * counts[name] / total : 0.0)
                }
                printf "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", \
                    pct("Promoter"), pct("Exon"), pct("Intron"), pct("Intergenic"), \
                    pct("Downstream"), pct("UTR5"), pct("UTR3"), pct("TTS"), pct("Other")
            }
        ' "$annot"
    )

    echo -e "${sample}\t${row}" >> "$HOMER_SUMMARY"
done

echo "HOMER annotation summary written to: $HOMER_SUMMARY"
# --- end HOMER summary ---


###############################################################################
# 3) Write MultiQC config (module order) + run MultiQC
###############################################################################
if [[ "$RUN_MULTIQC" -eq 1 ]]; then
  echo "=== Writing MultiQC config with custom module order & sample cleanup ==="
  MQC_CONFIG="${BASE}/multiqc_config.yaml"

  cat > "$MQC_CONFIG" <<'EOF'
# ---- Order of sections in the report ----
report_section_order:
  fastqc:
    order: 1000          # top
  cutadapt:
    order: 900
  bowtie2:
    order: 800
  samtools:
    order: 700
  preseq:
    order: 600
  macs2:
    order: 500
  phantompeakqualtools:
    order: 400
  homer:
    order: 300
  deeptools:
    order: 200           # still above misc stuff

# ---- Clean up sample names so each bio sample has ONE row ----
# Turn on regex mode:
sample_names_replace_regex: true

# Patterns to strip from sample names
# (Note double backslashes where needed if you ever use groups.)
sample_names_replace:
  "_R[12](_val_[12])?$": ""   # _R1, _R2, _R1_val_1, _R2_val_2 -> ""
  "_bowtie2$": ""             # alignment log suffix -> ""
  "_vs_.*$": ""               # comparison suffix like _vs_IgG -> ""

# ---- Ignore junk files so they don't show as samples ----
fn_ignore_files:
  - "*.sh"
  - "*multiqc_config.yaml"

# ---- Custom modules / plots ----
custom_data:
  frip_scores:
    file_format: "tsv"
    section_name: "FRiP scores (MACS2)"
    description: "Fraction of mapped reads falling in MACS2 peak regions for each sample."
    plot_type: "bargraph"
    pconfig:
      id: "frip_scores_plot"
      title: "FRiP scores"
      ylab: "FRiP"
      ymin: 0
      ymax: 1

  homer_peak_annotation:
    file_format: "tsv"
    section_name: "HOMER peak annotation"
    description: "Genomic distribution of MACS2 peaks annotated by HOMER annotatePeaks.pl."
    plot_type: "bargraph"
    pconfig:
      id: "homer_peak_annotation"
      title: "Peak annotation (HOMER)"
      ylab: "Percentage of peaks"
      xlab: "Sample"

sp:
  frip_scores:
    # Adjust path pattern to where frip_scores.tsv lives
    fn: "frip_scores.tsv"

  homer_peak_annotation:
    # Adjust if you put this file in a subdir
    fn: "homer_peak_annotation.tsv"
EOF

  echo "=== Running MultiQC over BASE: $BASE ==="
  multiqc "$BASE" -o "${BASE}/multiqc_all" -c "$MQC_CONFIG" --force
  echo "MultiQC report written to: ${BASE}/multiqc_all"

fi