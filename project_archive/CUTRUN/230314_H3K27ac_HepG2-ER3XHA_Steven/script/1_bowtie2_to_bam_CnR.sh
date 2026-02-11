#!/bin/bash

# Description: This script performs Bowtie2 alignment on paired-end trimmed FASTQ files, make tagdirectories, converts SAM files to sorted BAM files, index, and coverage the track of the sorted BAM files. And then it performs multiBamSummary, plotPCA, plotCorrelation (Spearman), and plotFingerprint on sorted.bam files.
# Input file format: Paired-end trimmed FASTQ files with names ending in *_R1_val_1.fq.gz and *_R2_val_2.fq.gz
# Last modified time: 2023-05-23

# Set parameters and directories
HOME="/mnt/home/digan"
BASE="${HOME}/seq/CUTRUN/230314_H3K27ac_HepG2-ER3XHA_Steven"
INPUT_DIR="${BASE}/cleandata"
OUTPUT_DIR="${BASE}/align"
BOWTIE2="/mnt/share/apps/bio/bowtie2-2.4.5/bowtie2"
HOMER="/opt/apps/bio/homer-4.11.1/bin"
SAMTOOLS="/opt/apps/bio/samtools-1.15.1/bin/samtools"
DEEPTOOLS="/opt/apps/bio/deeptools-3.5.1/bin"

# Set parameters for bowtie2
CPU=32
GENOME="hg19"
GENOME_INDEX="/mnt/share/archive/ref/align/bowtie2/${GENOME}/${GENOME}"

# Set parameters for bamCoverage
BINSIZE=10
GENOMESIZE=2864785220
NORMALIZE="RPGC"
IGNORE="chrM"

# Set regions and labels for bamqc plot
CHROMOSOME=19
LABELSTART="CnR-"
LABELEND="_23"

# Set the LD_LIBRARY_PATH to include the path to the libhts.so.3 file
export LD_LIBRARY_PATH="/mnt/share/apps/local/lib:$LD_LIBRARY_PATH"

# Make output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/bamqc
mkdir -p $OUTPUT_DIR/track

# Set the log file name with the output directory path
LOGFILE="$OUTPUT_DIR/1_bowtie2_to_bam_qc_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

# Echo the samples being processed
input_files=$(ls ${INPUT_DIR}/*.fq.gz)
echo "Samples being processed: $input_files"

# Loop through R1_val_1.fq.gz files in the input directory
for R1 in "$INPUT_DIR"/*_R1_val_1.fq.gz; do
  # Replace R1_val_1 with R2_val_2 in the file name to get the corresponding R2_val_2.fq.gz file
  R2="${R1%_R1_val_1.fq.gz}_R2_val_2.fq.gz"

  # Get the sample name from the R1 file name
  sample_name=$(basename "$R1" _R1_val_1.fq.gz)

  # Define output file names for SAM, BAM, and log
  output_sam="${OUTPUT_DIR}/${sample_name}.sam"
  output_tag="${OUTPUT_DIR}/${sample_name}_tag"
  output_bam="${OUTPUT_DIR}/${sample_name}.bam"
  output_sorted_bam="${OUTPUT_DIR}/${sample_name}_sorted.bam"
  output_bw="${OUTPUT_DIR}/track/${sample_name}.bw"
  log_file="${OUTPUT_DIR}/${sample_name}_bowtie2.log"

  # Run Bowtie2 with the specified parameters
  echo "Running Bowtie2 alignment for sample: $sample_name"
  nice $BOWTIE2 -x $GENOME_INDEX -1 "$R1" -2 "$R2" --very-sensitive -p $CPU -S "$output_sam" 2>"$log_file"

  # MakeTagDirectory for each sample using HOMER
  nice $HOMER/makeTagDirectory $output_tag -format sam -unique -tbp 1 $output_sam

  # Convert SAM to BAM, sort, and index BAM files
  echo "Running Samtools operations for sample: $sample_name"
  nice $SAMTOOLS view -bS "$output_sam" >"$output_bam"
  nice $SAMTOOLS sort "$output_bam" -o "$output_sorted_bam"
  nice $SAMTOOLS index "$output_sorted_bam"

  # Generate coverage track
  echo "Genearating coverage track for sample: $sample_name"
  nice $DEEPTOOLS/bamCoverage -b "$output_sorted_bam" -o "$output_bw" -bs $BINSIZE -v \
    --effectiveGenomeSize $GENOMESIZE --normalizeUsing $NORMALIZE --ignoreForNormalization $IGNORE

  # Check if SAM, BAM, and sorted BAM files exist, then delete SAM and unsorted BAM files, or print an error message
  if [[ -f $output_sam && -f $output_bam && -f $output_sorted_bam ]]; then
    echo "Deleting SAM and unsorted BAM files for sample: $sample_name"
    rm "$output_sam"
    rm "$output_bam"
  else
    echo "Error: One or more of the required files (SAM, BAM, or sorted BAM) not found for sample: $sample_name"
    break
  fi

  echo "Aglinmen and processing complete for sample: $sample_name"
done

echo "Alignments and processing complete! Start running BAM quality control"

# Generate label legends for the samples
labels=""
for f in "${OUTPUT_DIR}"/*_sorted.bam; do
  sample_name=$(basename "$f" _sorted.bam)
  label="${sample_name#*"$LABELSTART"}"
  label="${label%"$LABELEND"*}"
  labels+="$label "
done

# Run multiBamSummary
echo "Running multiBamSummary..."
nice ${DEEPTOOLS}/multiBamSummary bins -b ${OUTPUT_DIR}/*_sorted.bam --region $CHROMOSOME \
  --outFileName ${OUTPUT_DIR}/bamqc/multiBamSummary.npz

# Run plotPCA
echo "Running plotPCA..."
nice ${DEEPTOOLS}/plotPCA -in ${OUTPUT_DIR}/bamqc/multiBamSummary.npz --labels $labels \
  -o ${OUTPUT_DIR}/bamqc/PCA_plot.png

# Run plotCorrelation
echo "Running plotCorrelation (Spearman)..."
nice ${DEEPTOOLS}/plotCorrelation -in ${OUTPUT_DIR}/bamqc/multiBamSummary.npz --corMethod spearman \
  --skipZeros --plotTitle "Spearman Correlation" --whatToPlot heatmap \
  --colorMap RdYlBu_r --plotNumbers --labels $labels \
  -o ${OUTPUT_DIR}/bamqc/correlation_heatmap.png

# Run plotFingerprint
echo "Running plotFingerprint..."
nice ${DEEPTOOLS}/plotFingerprint -b ${OUTPUT_DIR}/*_sorted.bam --labels $labels \
  --region $CHROMOSOME -o ${OUTPUT_DIR}/bamqc/fingerprint_plot.png

echo "All analyses completed."