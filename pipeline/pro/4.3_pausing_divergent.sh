#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# 4.3_pausing_divergent.sh — PRO-seq pausing index + divergent eRNA pairs
#
# Why this step exists
#   PRO-seq's canonical quantitations are not ChIP peak calls:
#     (a) Pausing index (PI) = promoter tag density / gene body tag density.
#         A high PI flags Pol II paused at the promoter; changes under
#         treatment capture pause release (e.g. E2 induction → PI drops at
#         induced genes).
#     (b) Divergent transcription: nearly every active promoter initiates
#         bidirectionally, producing a short antisense unstable RNA upstream
#         of the sense transcript. Divergent-pair detection flags active
#         promoters and is a useful enhancer/promoter classifier.
#
# Inputs
#   - ${BAM_DIR}/*_sorted.bam  (with .bai)
#   - ${GENE_ANNOT_BED}        BED(12 ok — first 6 cols used: chr start end name score strand)
#   - ${HOMER_PEAK_DIR}/*.annotatePeaks.txt — optional; used for divergent
#     detection (findPeaks -style groseq output from step 4.2). If missing,
#     divergent detection is skipped for that sample.
#
# Outputs (in ${PAUSING_DIR}/)
#   promoter.bed, genebody.bed             one-time cached feature BEDs
#   ${sample}_promoter_counts.bed           bedtools multicov counts
#   ${sample}_genebody_counts.bed
#   ${sample}_pausing_index.tsv             gene_id, strand, prom_density, gb_density, PI
#   ${sample}_divergent_pairs.tsv           sense_id, antisense_id, distance_bp
#   pausing_index_summary_mqc.tsv           median/quartile PI per sample for MultiQC
#
# Strand handling
#   PRO-seq R1 aligns to the OPPOSITE strand of the nascent RNA. For BAM-based
#   strand-aware counting we want reads whose aligned strand is OPPOSITE to the
#   feature (gene) strand. bedtools multicov -S does this; multicov -s would
#   count sense reads. Respect PROSEQ_FLIP_STRAND — if 0, use -s (read strand
#   == RNA strand, i.e. pre-flipped chemistry or unflipped by user choice).
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/0_config.sh"

log_start "4.3_pausing_divergent"

check_cmd bedtools
check_cmd awk
check_cmd samtools

echo "=== STEP 4.3: Pausing index + divergent eRNA ==="
echo "[INFO] BAM_DIR:            $BAM_DIR"
echo "[INFO] PAUSING_DIR:        $PAUSING_DIR"
echo "[INFO] GENE_ANNOT_BED:     $GENE_ANNOT_BED"
echo "[INFO] TSS_WINDOW:         $TSS_WINDOW (bp around TSS)"
echo "[INFO] GB_START:           $GB_START (bp downstream of TSS for gene body)"
echo "[INFO] GB_MIN_LEN:         $GB_MIN_LEN (minimum gene length)"
echo "[INFO] DIVERGENT_WINDOW:   $DIVERGENT_WINDOW (bp for antisense pair search)"
echo "[INFO] PROSEQ_FLIP_STRAND: ${PROSEQ_FLIP_STRAND:-0}"

if [[ ! -f "$GENE_ANNOT_BED" ]]; then
  echo "[ERROR] GENE_ANNOT_BED not found: $GENE_ANNOT_BED"
  exit 1
fi

mkdir -p "$PAUSING_DIR"

# Reads-vs-feature strand flag for bedtools multicov.
# PROSEQ_FLIP_STRAND=1 → count reads on OPPOSITE strand of gene (aligned strand
#                       = antisense of RNA, and we want "nascent RNA" reads).
# PROSEQ_FLIP_STRAND=0 → count reads on SAME strand as gene.
if [[ "${PROSEQ_FLIP_STRAND:-0}" -eq 1 ]]; then
  COUNT_STRAND_FLAG="-S"
else
  COUNT_STRAND_FLAG="-s"
fi
echo "[INFO] bedtools multicov strand flag: $COUNT_STRAND_FLAG"

###############################################################################
# 1. Build promoter / gene body feature BEDs from the annotation
###############################################################################

PROM_BED="${PAUSING_DIR}/promoter.bed"
GB_BED="${PAUSING_DIR}/genebody.bed"

if [[ ! -s "$PROM_BED" || ! -s "$GB_BED" ]]; then
  echo "[STEP4.3] Building promoter.bed (TSS ± ${TSS_WINDOW}) and genebody.bed (TSS+${GB_START}..TTS)"
  # Use only first 6 columns; keep one interval per transcript ID.
  awk -v W="$TSS_WINDOW" -v OFS='\t' '
    NF >= 6 {
      name = $4; score = $5; strand = $6
      if (strand == "+") { s = $2 - W; if (s < 0) s = 0; e = $2 + W; print $1, s, e, name, score, strand }
      else if (strand == "-") { s = $3 - W; if (s < 0) s = 0; e = $3 + W; print $1, s, e, name, score, strand }
    }
  ' "$GENE_ANNOT_BED" | sort -k1,1 -k2,2n > "$PROM_BED"

  awk -v GB="$GB_START" -v MIN="$GB_MIN_LEN" -v OFS='\t' '
    NF >= 6 && ($3 - $2) >= MIN {
      name = $4; score = $5; strand = $6
      if (strand == "+") { print $1, $2 + GB, $3, name, score, strand }
      else if (strand == "-") { print $1, $2, $3 - GB, name, score, strand }
    }
  ' "$GENE_ANNOT_BED" | sort -k1,1 -k2,2n > "$GB_BED"

  echo "[STEP4.3] promoter.bed: $(wc -l < "$PROM_BED") intervals"
  echo "[STEP4.3] genebody.bed: $(wc -l < "$GB_BED") intervals"
fi

###############################################################################
# 2. Per-sample pausing index
###############################################################################

# Split summaries so MultiQC custom_content renders gene counts and PI quartiles
# in their own tables — combining them put N_genes (~10^5) on the same axis as
# Median/Q1/Q3 PI (~10^0) and crushed the PI bars to invisibility.
PI_QUART_SUMMARY="${PAUSING_DIR}/pausing_index_quartiles_mqc.tsv"
PI_NGENES_SUMMARY="${PAUSING_DIR}/pausing_index_ngenes_mqc.tsv"
echo -e "Sample\tMedian_PI\tQ1_PI\tQ3_PI" > "$PI_QUART_SUMMARY"
echo -e "Sample\tN_genes" > "$PI_NGENES_SUMMARY"
# Remove the legacy combined file so MultiQC's _mqc.tsv auto-detection
# doesn't surface a stale section alongside the new split tables.
rm -f "${PAUSING_DIR}/pausing_index_summary_mqc.tsv"

sorted_bams=( "${BAM_DIR}"/*_sorted.bam )
if (( ${#sorted_bams[@]} == 0 )); then
  echo "[WARN] No *_sorted.bam in $BAM_DIR — nothing to do."
  exit 0
fi

for bam in "${sorted_bams[@]}"; do
  sample="$(basename "$bam" _sorted.bam)"
  prom_cnt="${PAUSING_DIR}/${sample}_promoter_counts.bed"
  gb_cnt="${PAUSING_DIR}/${sample}_genebody_counts.bed"
  pi_tsv="${PAUSING_DIR}/${sample}_pausing_index.tsv"

  echo "[STEP4.3] $sample: counting reads in promoter and gene body"

  if [[ ! -s "$prom_cnt" ]]; then
    bedtools multicov -bams "$bam" -bed "$PROM_BED" $COUNT_STRAND_FLAG > "$prom_cnt"
  fi
  if [[ ! -s "$gb_cnt" ]]; then
    bedtools multicov -bams "$bam" -bed "$GB_BED" $COUNT_STRAND_FLAG > "$gb_cnt"
  fi

  # Join on transcript ID (column 4) and compute PI.
  # prom_cnt/gb_cnt columns: chr start end name score strand count
  # density = count / (end - start); PI = (prom_density + eps) / (gb_density + eps); eps avoids /0
  echo -e "gene_id\tstrand\tprom_len\tprom_count\tprom_density\tgb_len\tgb_count\tgb_density\tpausing_index" > "$pi_tsv"

  awk -v OFS='\t' '{ print $4, $3-$2, $7, $6 }' "$prom_cnt" | sort -k1,1 > "${pi_tsv}.prom.tmp"
  awk -v OFS='\t' '{ print $4, $3-$2, $7, $6 }' "$gb_cnt"   | sort -k1,1 > "${pi_tsv}.gb.tmp"

  # join -t tab; only genes present in both (i.e. gene_length >= GB_MIN_LEN).
  # After join -1 1 -2 1 on two 4-col files (id, length, count, strand) the
  # merged line is: $1=id, $2=prom_len, $3=prom_count, $4=prom_strand,
  #                 $5=gb_len, $6=gb_count, $7=gb_strand
  join -t $'\t' -1 1 -2 1 "${pi_tsv}.prom.tmp" "${pi_tsv}.gb.tmp" \
    | awk -v OFS='\t' -v EPS=1e-9 '
      {
        gene_id=$1; prom_len=$2; prom_count=$3; strand=$4;
        gb_len=$5; gb_count=$6;
        if (prom_len <= 0 || gb_len <= 0) next
        prom_d = prom_count / prom_len
        gb_d   = gb_count   / gb_len
        pi = (prom_d + EPS) / (gb_d + EPS)
        printf "%s\t%s\t%d\t%d\t%.6g\t%d\t%d\t%.6g\t%.4f\n", gene_id, strand, prom_len, prom_count, prom_d, gb_len, gb_count, gb_d, pi
      }' >> "$pi_tsv"
  rm -f "${pi_tsv}.prom.tmp" "${pi_tsv}.gb.tmp"

  # Summary row: restrict to genes with reads in BOTH promoter AND gene body.
  # Without this, genes with zero gene-body reads blow up PI to ≈ prom_d/EPS ≈ 1e9
  # and poison the median. The per-gene TSV keeps every gene (unfiltered) so
  # downstream analyses are free to apply their own thresholds.
  read -r n_genes med_pi q1_pi q3_pi < <(
    tail -n +2 "$pi_tsv" | awk '$4 > 0 && $7 > 0 {print $9}' | sort -g | awk '
      function q(p,   x, i, f) {
        x = p * (n - 1) + 1
        i = int(x); f = x - i
        if (i < 1) return a[1]
        if (i >= n) return a[n]
        return a[i] + f*(a[i+1]-a[i])
      }
      { a[NR]=$1 }
      END {
        n=NR
        if (n == 0) { print 0, "NA", "NA", "NA"; exit }
        printf "%d %.4f %.4f %.4f\n", n, q(0.5), q(0.25), q(0.75)
      }'
  )
  echo -e "${sample}\t${med_pi}\t${q1_pi}\t${q3_pi}" >> "$PI_QUART_SUMMARY"
  echo -e "${sample}\t${n_genes}" >> "$PI_NGENES_SUMMARY"
  echo "  → genes=${n_genes}  median PI=${med_pi}  Q1/Q3=${q1_pi}/${q3_pi}"
done

###############################################################################
# 3. Divergent eRNA pairs — from step 4.2 annotatePeaks output (groseq style)
###############################################################################

echo "[STEP4.3] Divergent eRNA pair detection (DIVERGENT_WINDOW=${DIVERGENT_WINDOW} bp)"

extract_5p_bed_from_annotate() {
  # annotatePeaks.txt: col1=PeakID, col2=Chr, col3=Start, col4=End, col5=Strand
  # 5′ end per strand: + → Start, - → End
  local ann="$1"
  tail -n +2 "$ann" \
    | awk -v OFS='\t' 'NF>=5 && $2 != "" {
        if ($5=="+") print $2, $3, $3+1, $1, ".", $5
        else if ($5=="-") print $2, ($4>0 ? $4-1 : 0), $4, $1, ".", $5
      }'
}

for ann in "${HOMER_PEAK_DIR}"/*.annotatePeaks.txt; do
  [[ -f "$ann" ]] || continue
  name="$(basename "$ann" .annotatePeaks.txt)"
  out="${PAUSING_DIR}/${name}_divergent_pairs.tsv"
  tmp_5p="${PAUSING_DIR}/${name}_5p.bed"

  echo "  · $name"
  extract_5p_bed_from_annotate "$ann" | sort -k1,1 -k2,2n > "$tmp_5p"
  plus_5p="${tmp_5p%.bed}.plus.bed"
  minus_5p="${tmp_5p%.bed}.minus.bed"
  awk '$6=="+"' "$tmp_5p" > "$plus_5p"
  awk '$6=="-"' "$tmp_5p" > "$minus_5p"

  # A + strand transcript at 5′=p has a divergent pair if a - strand transcript's 5′=m
  # lies within DIVERGENT_WINDOW upstream (m ≤ p, p - m ≤ W). Emit the closest such pair.
  echo -e "sense_id\tsense_chr\tsense_5p\tantisense_id\tantisense_chr\tantisense_5p\tdistance_bp" > "$out"
  if [[ -s "$plus_5p" && -s "$minus_5p" ]]; then
    bedtools window -a "$plus_5p" -b "$minus_5p" -l "$DIVERGENT_WINDOW" -r 0 \
      | awk -v OFS='\t' '
        # +strand a cols 1-6, -strand b cols 7-12; antisense 5p is at $8 (start of single-bp feature)
        { d = $2 - $8; if (d >= 0) print $4, $1, $2, $10, $7, $8, d }' \
      | sort -k1,1 -k7,7n \
      | awk -v OFS='\t' '!seen[$1]++' \
      >> "$out"
  fi

  rm -f "$tmp_5p" "$plus_5p" "$minus_5p"
  n_div=$(tail -n +2 "$out" | wc -l | tr -d ' ')
  echo "    divergent pairs: $n_div"
done

echo "=== STEP 4.3 (pausing + divergent) complete ==="
