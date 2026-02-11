#!/bin/bash

# Description: This script trims the input FASTQ files using Trim Galore, then performs FastQC quality control and MultiQC summary report generation.
# Input file format: Paired-end FASTQ files with names ending in *_R1.fastq.gz and *_R2.fastq.gz (e.g., sample1_R1.fastq.gz and sample1_R2.fastq.gz)
# Last modified time: 2023-07-27

# Set parameters and directories
HOME="/mnt/home/digan"
BASE="${HOME}/seq/CUTRUN/230530_CUTRUN_NCAPG-H-H2_HepG2_ER3XHA_Steven"
TRIM_GALORE="/opt/apps/bio/TrimGalore-0.6.10/trim_galore"
CUTADAPT="/opt/apps/bio/cutadapt-4.1/bin/cutadapt"
FASTQC="/opt/apps/bio/fastqc-0.11.9/fastqc"
MULTIQC="/opt/apps/bio/multiqc-1.12/bin/multiqc"
INPUT_DIR="${BASE}/data"
OUTPUT_DIR="${BASE}/cleandata"

# Set parameters for trimming program
STRINGENCY=5
CPU=8

# Make output directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/fastqc
mkdir -p ${OUTPUT_DIR}/fastqc/multiqc

# Set the log file name with the output directory path
LOGFILE="$OUTPUT_DIR/0_trim_galore_qc_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

# Echo the samples being processed
input_files=$(ls ${INPUT_DIR}/*.fastq.gz)
echo "Samples being processed: $input_files"

# Loop through R1.fastq.gz files in the input directory
for R1 in "$INPUT_DIR"/*R1.fastq.gz; do
  # Replace R1 with R2 in the file name to get the corresponding R2.fastq.gz file
  R2="${R1%R1.fastq.gz}R2.fastq.gz"

  # Get the sample name from the R1 file name
  sample_name=$(basename "$R1" _R1.fastq.gz)

  echo "Running Galore for sample: $sample_name"
  # Run Trim Galore with the specified parameters
  nice $TRIM_GALORE --paired --output_dir $OUTPUT_DIR -q 25 --phred33 \
    --stringency $STRINGENCY --length 36 --path_to_cutadapt $CUTADAPT --gzip "$R1" "$R2"
done

echo "Running FastQC for all samples"
# Run FastQC with the specified parameters
nice $FASTQC --outdir ${OUTPUT_DIR}/fastqc --threads $CPU ${OUTPUT_DIR}/*_val_*.fq.gz

echo "Running MultiQC for all samples"
# Run MultiQC
nice $MULTIQC ${OUTPUT_DIR}/fastqc -o ${OUTPUT_DIR}/fastqc/multiqc

echo "Analyses complete!"
