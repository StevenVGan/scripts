#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# CUT&RUN UPSTREAM PIPELINE
# (trimming → QC → alignment → BAM processing → HOMER tags → MACS3 peaks + FRiP)
#
# - Assumes conda env "bio" is active and tools are on PATH
# - Input format:
#     data/*.fastq.gz named like: sample_R1.fastq.gz / sample_R2.fastq.gz
# - Trim output:
#     cleandata/*_R1_val_1.fq.gz / *_R2_val_2.fq.gz
# - Alignment output:
#     align/*_sorted.bam + bigWig tracks + optional deepTools multiBamSummary
# - MACS3 groups:
#     defined in macs3_groups.tsv (see format in CONFIG)
# - This script:
#     * NO MultiQC (handled in a separate QC script)
#     * Deletes trimmed FASTQs after alignment
#     * Builds HOMER tag directories from sorted BAMs
#     * Calls MACS3 (narrow vs broad based on TF/Histone)
#     * Computes FRiP per MACS3 group into peaks/frip_scores.tsv
###############################################################################

############################ CONFIG SECTION ###################################

# TODO: This code is temporarly using macs3 via 'conda run -n macs3 macs3 ...' instead of macs3 directly on PATH. This is because of a zerodivision bug of macs3. Need to fix this properly later. 

# TODO: I have update to macs3, no need to keep the old macs3 conda env. Later we can just use macs3 directly on PATH.

# Which steps to run (1 = yes, 0 = no)
RUN_TRIMMING=0
RUN_ALIGNMENT=1
RUN_MACS3=1

# Project base directory
BASE="${HOME}/work/seq/CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka"

# Directories (relative to BASE)
RAW_DIR="${BASE}/data"           # symlinked + renamed FASTQs
TRIM_DIR="${BASE}/cleandata"
ALIGN_DIR="${BASE}/align"
TRACK_DIR="${ALIGN_DIR}/track"
BAM_DIR="${ALIGN_DIR}/bam"
BAMQC_DIR="${ALIGN_DIR}/bamqc"
TAG_DIR="${ALIGN_DIR}/tags"
PEAK_DIR="${BASE}/peaks"
MACS3_DIR="${PEAK_DIR}/macs3"
LOG_DIR="${BASE}/logs"

# MACS3 group definition file
# Expected format per line (tab- or space-separated):
#
#   ip_bam   control_bam   name   type
#
# - ip_bam & control_bam can be:
#     * just a filename (resolved relative to BAM_DIR)
#     * OR a full/relative path anywhere on the system
# - control_bam:
#     * "-" or "none"/"NA" → no control
# - type:
#     * "TF" or "Histone" (case-insensitive)
#     * if missing, defaults to "TF"
#
# Examples:
#   SG-339_CnR_ERa_MCF7_Veh_E2_rep1_sorted.bam SG-348_CnR_IgG_MCF7_Veh_ICI_sorted.bam ERa_E2_rep1_vs_IgG TF
#   SG-347_CnR_H3K4me3_MCF7_Veh_E2_sorted.bam  -                                  H3K4me3_E2_noControl Histone
MACS3_GROUPS_FILE="${BASE}/macs3_groups.tsv"

# Trimming parameters
TRIM_STRINGENCY=5
TRIM_MIN_LENGTH=36
TRIM_QUAL=25
TRIM_CPU=8

# Bowtie2 / alignment parameters
BT2_CPU=12
GENOME="hg38"
GENOME_INDEX="/mnt/share/archive/bkup/ref/align/bowtie2/${GENOME}/${GENOME}"

# BAM coverage / deeptools parameters (effective genome size for hg38)
BINSIZE=10
GENOMESIZE=2913022398     # hg38 effective genome size (approx)
NORMALIZE="RPGC"
IGNORE_CHR="chrM"

# deepTools multiBamSummary (NPZ only; plots handled later)
BAMQC_REGION="19"   # region for multiBamSummary if desired

# MACS3 peak calling parameters
MACS3_GENOMESIZE="hs"
MACS3_FORMAT="BAMPE"
MACS3_FDR=0.05
MACS3_BROAD_CUTOFF=0.10            # used for Histone samples when calling broad peaks
BLACKLIST="${HOME}/work/ref/blacklist/${GENOME}/${GENOME}-blacklist.v2.bed"

######################## END CONFIG SECTION ###################################

# Helper: check that required commands exist
check_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command '$cmd' not found in PATH. Activate the correct conda env (e.g. 'bio')." >&2
    exit 1
  fi
}

echo "=== CUT&RUN upstream pipeline starting ==="

# Check tools (no explicit paths; rely on conda 'bio' env)
for c in trim_galore cutadapt fastqc bowtie2 samtools bamCoverage \
         multiBamSummary bedtools makeTagDirectory; do
  check_cmd "$c"
done

# Create directories
mkdir -p "$TRIM_DIR" "$ALIGN_DIR" "$BAM_DIR" "$TRACK_DIR" "$BAMQC_DIR" "$TAG_DIR" "$PEAK_DIR" "$MACS3_DIR" "$LOG_DIR"

# Global log file
PIPELINE_LOG="${LOG_DIR}/1_trim_bowtie2_macs3_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to: $PIPELINE_LOG"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

echo "[INFO] BASE: $BASE"
echo "[INFO] RAW_DIR: $RAW_DIR"
echo "[INFO] TRIM_DIR: $TRIM_DIR"
echo "[INFO] ALIGN_DIR: $ALIGN_DIR"
echo "[INFO] PEAK_DIR: $PEAK_DIR"

###############################################################################
# STEP 1: Trimming (Trim Galore) + FastQC   (NO MultiQC here)
###############################################################################
if [[ "$RUN_TRIMMING" -eq 1 ]]; then
  echo "=== STEP 1: Trimming + FastQC ==="
  echo "[INFO] Looking for raw FASTQ files in: $RAW_DIR"

  ls "$RAW_DIR"/*.fastq.gz 2>/dev/null || echo "[WARN] No FASTQ files found in $RAW_DIR"

  # Trim Galore
  for R1 in "$RAW_DIR"/*_R1.fastq.gz; do
    [[ -e "$R1" ]] || { echo "[WARN] No *_R1.fastq.gz files found in $RAW_DIR"; break; }

    R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    if [[ ! -e "$R2" ]]; then
      echo "[WARN] Paired file not found for $R1 (expected $R2), skipping."
      continue
    fi

    sample_name=$(basename "$R1" _R1.fastq.gz)
    echo "[STEP1] Trimming sample: $sample_name"

    trim_galore \
      --paired \
      --output_dir "$TRIM_DIR" \
      -q "$TRIM_QUAL" \
      --phred33 \
      --stringency "$TRIM_STRINGENCY" \
      --length "$TRIM_MIN_LENGTH" \
      --gzip \
      "$R1" "$R2"
  done

  echo "[STEP1] Running FastQC on trimmed reads..."
  mkdir -p "$TRIM_DIR/fastqc"

  fastqc \
    --outdir "$TRIM_DIR/fastqc" \
    --threads "$TRIM_CPU" \
    "$TRIM_DIR"/*_val_*.fq.gz

  echo "=== STEP 1 complete ==="
else
  echo "=== STEP 1 skipped (RUN_TRIMMING=0) ==="
fi

###############################################################################
# STEP 2: Bowtie2 → BAM + coverage + HOMER tag dirs + multiBamSummary NPZ
###############################################################################
if [[ "$RUN_ALIGNMENT" -eq 1 ]]; then
  echo "=== STEP 2: Alignment + BAM processing + HOMER tag directories ==="
  echo "[INFO] Using Bowtie2 index: $GENOME_INDEX"

  ls "$TRIM_DIR"/*_R1_val_1.fq.gz 2>/dev/null || echo "[WARN] No trimmed FASTQ files found in $TRIM_DIR"

  for R1 in "$TRIM_DIR"/*_R1_val_1.fq.gz; do
    [[ -e "$R1" ]] || { echo "[WARN] No *_R1_val_1.fq.gz files found in $TRIM_DIR"; break; }

    R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"
    if [[ ! -e "$R2" ]]; then
      echo "[WARN] Paired trimmed file not found for $R1 (expected $R2), skipping."
      continue
    fi

    sample_name=$(basename "$R1" _R1_val_1.fq.gz)
    echo "[STEP2] Processing sample: $sample_name"

    sam_file="${ALIGN_DIR}/${sample_name}.sam"
    bam_file="${BAM_DIR}/${sample_name}.bam"
    sorted_bam="${BAM_DIR}/${sample_name}_sorted.bam"
    
    bw_file="${TRACK_DIR}/${sample_name}.bw"
    bt2_log="${BAM_DIR}/${sample_name}_bowtie2.log"

    echo "[STEP2] Bowtie2 alignment..."
    bowtie2 \
      -x "$GENOME_INDEX" \
      -1 "$R1" \
      -2 "$R2" \
      --very-sensitive \
      -p "$BT2_CPU" \
      -S "$sam_file" \
      2> "$bt2_log"

    echo "[STEP2] Samtools: SAM → BAM → sorted BAM → index..."
    samtools view -bS "$sam_file" > "$bam_file"
    samtools sort "$bam_file" -o "$sorted_bam"
    samtools index "$sorted_bam"
    # Simple per-contig mapping stats for MultiQC (Mapped reads per contig)
    samtools idxstats "$sorted_bam" > "${BAMQC_DIR}/${sample_name}.idxstats.txt"

    echo "[STEP2] bamCoverage: generating bigWig track..."
    bamCoverage \
      -b "$sorted_bam" \
      -o "$bw_file" \
      -bs "$BINSIZE" \
      --effectiveGenomeSize "$GENOMESIZE" \
      --normalizeUsing "$NORMALIZE" \
      --ignoreForNormalization "$IGNORE_CHR" \
      --ignoreDuplicates

    echo "[STEP2] HOMER makeTagDirectory for $sample_name..."
    makeTagDirectory "${TAG_DIR}/${sample_name}" "$sam_file" -format sam -tbp 1

    if [[ -f "$sam_file" && -f "$bam_file" && -f "$sorted_bam" ]]; then
      echo "[STEP2] Cleaning up SAM and unsorted BAM for $sample_name"
      rm -f "$sam_file" "$bam_file"
    else
      echo "[ERROR] Missing SAM/BAM/sorted BAM for $sample_name; not deleting intermediates."
    fi

    echo "[STEP2] Alignment + processing complete for sample: $sample_name"
  done

  echo "[STEP2] deepTools multiBamSummary on all sorted BAMs (NPZ only)..."
  sorted_bams=( "$BAM_DIR"/*_sorted.bam )
  if (( ${#sorted_bams[@]} == 0 )); then
    echo "[WARN] No *_sorted.bam files found in $BAM_DIR; skipping multiBamSummary."
  else
    # Build labels as simple sample names (basename without _sorted.bam)
    labels=""
    for f in "${sorted_bams[@]}"; do
      s=$(basename "$f" _sorted.bam)
      labels+="$s "
    done

    QC_NPZ="${BAMQC_DIR}/multiBamSummary.npz"

    echo "[STEP2] multiBamSummary..."
    multiBamSummary bins \
      -b "${sorted_bams[@]}" \
      --region "$BAMQC_REGION" \
      --outFileName "$QC_NPZ"
    # No plotPCA / plotCorrelation / plotFingerprint here;
    # they will be handled later in a dedicated QC/MultiQC script.
  fi

  echo "[STEP2] Deleting trimmed FASTQs to save space..."
  rm -f "$TRIM_DIR"/*_val_*.fq.gz || true

  echo "=== STEP 2 complete ==="
else
  echo "=== STEP 2 skipped (RUN_ALIGNMENT=0) ==="
fi

###############################################################################
# STEP 3: MACS3 peaks (flexible controls, TF vs Histone, FRiP per group)
###############################################################################
if [[ "$RUN_MACS3" -eq 1 ]]; then
  echo "=== STEP 3: MACS3 peak calling (no bamCompare) ==="
  echo "[INFO] MACS3 groups file: $MACS3_GROUPS_FILE"

  if [[ ! -f "$MACS3_GROUPS_FILE" ]]; then
    echo "[WARN] MACS3 groups file not found; skipping MACS3 step."
  else
    # FRiP results for all groups
    FRIP_TSV="${MACS3_DIR}/frip_scores.tsv"
    echo -e "sample\treads_in_peaks\ttotal_mapped\tfrip" > "$FRIP_TSV"

    while read -r ip_spec ctrl_spec name type _rest; do
      # Skip blank lines or comment lines
      [[ -z "${ip_spec:-}" ]] && continue
      [[ "$ip_spec" =~ ^# ]] && continue

      # Default name from ip if missing
      if [[ -z "${name:-}" ]]; then
        name="$(basename "$ip_spec" .bam)"
      fi

      # Default type to TF if missing
      if [[ -z "${type:-}" ]]; then
        type="TF"
      fi

      # Normalise type case
      type_upper="$(echo "$type" | tr '[:lower:]' '[:upper:]')"

      # Resolve IP BAM path:
      #   - If ip_spec exists as given, use it
      #   - Else, try BAM_DIR/ip_spec
      if [[ -f "$ip_spec" ]]; then
        ip_file="$ip_spec"
      elif [[ -f "${BAM_DIR}/${ip_spec}" ]]; then
        ip_file="${BAM_DIR}/${ip_spec}"
      else
        echo "[WARN] IP BAM not found for spec '$ip_spec' (checked as-is and in $BAM_DIR); skipping group $name"
        continue
      fi

      # Determine whether we have a control
      use_control=1
      if [[ -z "${ctrl_spec:-}" ]] || [[ "$ctrl_spec" == "-" ]] || [[ "$ctrl_spec" =~ ^(none|NONE|NA)$ ]]; then
        use_control=0
      fi

      control_file=""
      if [[ "$use_control" -eq 1 ]]; then
        # Resolve control path similarly:
        if [[ -f "$ctrl_spec" ]]; then
          control_file="$ctrl_spec"
        elif [[ -f "${BAM_DIR}/${ctrl_spec}" ]]; then
          control_file="${BAM_DIR}/${ctrl_spec}"
        else
          echo "[WARN] Control BAM not found for spec '$ctrl_spec' (checked as-is and in $BAM_DIR); skipping group $name"
          continue
        fi
      fi

      echo "[STEP3] MACS3 group: NAME=$name  TYPE=$type_upper"
      echo "        IP:   $ip_file"
      if [[ "$use_control" -eq 1 ]]; then
        echo "        CTRL: $control_file"
      else
        echo "        CTRL: (none)"
      fi

      # Build macs3 command
      macs_cmd=( conda run -n macs3 macs3 callpeak
        -t "$ip_file"
        -g "$MACS3_GENOMESIZE"
        -q "$MACS3_FDR"
        --format "$MACS3_FORMAT"
        --outdir "$MACS3_DIR"
        --name "$name"
      )

      # Add control if present
      if [[ "$use_control" -eq 1 ]]; then
        macs_cmd+=( -c "$control_file" )
      fi

      # Use broad mode for Histone marks
      if [[ "$type_upper" == "HISTONE" ]]; then
        macs_cmd+=( --broad --broad-cutoff "$MACS3_BROAD_CUTOFF" )
      fi

      # Run MACS3 using macs env
      "${macs_cmd[@]}"

      narrowPeak_file="${MACS3_DIR}/${name}_peaks.narrowPeak"
      broadPeak_file="${MACS3_DIR}/${name}_peaks.broadPeak"
      filtered_bed="${MACS3_DIR}/${name}_filtered.bed"

      if [[ ! -f "$narrowPeak_file" && ! -f "$broadPeak_file" ]]; then
        echo "[WARN] No peak file (narrow or broad) found for $name, skipping blacklist filtering + FRiP."
        continue
      fi

      echo "[STEP3] Filtering peaks with blacklist: $BLACKLIST"
      if [[ -f "$narrowPeak_file" ]]; then
        bedtools intersect \
          -a "$narrowPeak_file" \
          -b "$BLACKLIST" \
          -v \
          > "$filtered_bed"
      elif [[ -f "$broadPeak_file" ]]; then
        bedtools intersect \
          -a "$broadPeak_file" \
          -b "$BLACKLIST" \
          -v \
          > "$filtered_bed"
      fi

      echo "[STEP3] MACS3 + filtering complete for group: $name"

      # ---------------- FRiP calculation (Fraction of reads in peaks) ----------------
      if command -v samtools >/dev/null 2>&1 && command -v bedtools >/dev/null 2>&1; then
        echo "[STEP3] Calculating FRiP for $name ..."

        # Total mapped reads (exclude unmapped [0x4], secondary [0x100], supplementary [0x800])
        total_mapped=$(samtools view -c -F 260 "$ip_file")

        # Reads that overlap peaks at least once
        reads_in_peaks=$(bedtools intersect -a "$ip_file" -b "$filtered_bed" -u \
          | samtools view -c -)

        # Compute FRiP as a floating-point ratio
        frip=$(awk -v inpeaks="$reads_in_peaks" -v total="$total_mapped" \
          'BEGIN { if (total > 0) printf "%.4f", inpeaks/total; else print "0.0000" }')

        echo -e "${name}\t${reads_in_peaks}\t${total_mapped}\t${frip}" >> "$FRIP_TSV"
      else
        echo "[WARN] samtools or bedtools not found, skipping FRiP for $name"
      fi
      # ---------------- end FRiP calculation ----------------

    done < "$MACS3_GROUPS_FILE"
  fi

  echo "=== STEP 3 complete ==="
else
  echo "=== STEP 3 skipped (RUN_MACS3=0) ==="
fi

echo "=== CUT&RUN upstream pipeline finished successfully ==="
