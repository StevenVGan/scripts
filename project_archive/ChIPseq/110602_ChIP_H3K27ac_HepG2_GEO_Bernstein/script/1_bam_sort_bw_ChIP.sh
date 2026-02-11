#!/bin/bash

# Description: This script performs Bowtie2 alignment on paired-end trimmed FASTQ files, converts SAM files to sorted BAM files, index, and coverage the track of the sorted BAM files. And then it performs multiBamSummary, plotPCA, plotCorrelation (Spearman), and plotFingerprint on sorted.bam files.
# Input file format: Paired-end trimmed FASTQ files with names ending in *_R1_val_1.fq.gz and *_R2_val_2.fq.gz
# Last modified time: 2023-05-23
# Modification: This modified script only process .bam files in the output directory, without running Bowtie2 alignment.

# Set parameters and directories
HOME="/mnt/home/digan"
BASE="${HOME}/seq/ChIPseq/110602_ChIP_H3K27ac_HepG2_GEO_Bernstein"
INPUT_DIR="${BASE}/align"
OUTPUT_DIR="${BASE}/align"
SAMTOOLS="/opt/apps/bio/samtools-1.15.1/bin/samtools"
DEEPTOOLS="/opt/apps/bio/deeptools-3.5.1/bin"

# Set parameters for bamCoverage
BINSIZE=10
GENOMESIZE=2864785220
NORMALIZE="RPGC"
IGNORE="chrM"

# Set regions and labels for bamqc plot
CHROMOSOME=19
LABELSTART="ChIP_"
LABELEND="_GEO"

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
input_files=$(ls ${INPUT_DIR}/*.bam)
echo "Samples being processed: $input_files"

# Loop through all the .bam files in the output directory
for BAM in "$OUTPUT_DIR"/*.bam; do

  # Get the sample name from the .bam file name
  sample_name=$(basename "$BAM" .bam)

  # Define output file names for SAM, BAM, and log
  output_sorted_bam="${OUTPUT_DIR}/${sample_name}_sorted.bam"
  output_bw="${OUTPUT_DIR}/track/${sample_name}.bw"
  
  # Convert SAM to BAM, sort, and index BAM files
  echo "Running Samtools operations for sample: $sample_name"
  nice $SAMTOOLS sort "$BAM" -o "$output_sorted_bam"
  nice $SAMTOOLS index "$output_sorted_bam"

  # Generate coverage track
  echo "Genearating coverage track for sample: $sample_name"
  nice $DEEPTOOLS/bamCoverage -b "$output_sorted_bam" -o "$output_bw" -bs $BINSIZE -v \
    --effectiveGenomeSize $GENOMESIZE --normalizeUsing $NORMALIZE --ignoreForNormalization $IGNORE

  echo "Processing complete for sample: $sample_name"
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