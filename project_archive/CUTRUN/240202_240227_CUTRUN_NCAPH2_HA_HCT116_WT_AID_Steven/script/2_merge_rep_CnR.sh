#!/bin/bash

# Description: This script merges sorted BAM files based on user input
# Input file format: Sorted BAM files with names ending in *_sorted.bam
# Last modified time: 2024-04-04

# Set parameters and directories
HOME="/mnt/home/digan/work"
BASE="${HOME}/seq/CUTRUN/240202_240227_CUTRUN_NCAPH2_HA_HCT116_WT_AID_Steven"
INPUT_DIR="${BASE}/align"
PICARD="/opt/apps/bio/picard-2.18.1/picard.jar"
SAMTOOLS="/opt/apps/bio/samtools-1.15.1/bin/samtools"

# Set the LD_LIBRARY_PATH to include the path to the libhts.so.3 file
export LD_LIBRARY_PATH="/mnt/share/apps/local/lib:$LD_LIBRARY_PATH"

# Set the log file name with the output directory path
LOGFILE="$INPUT_DIR/2_merge_rep_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

# Print all sorted BAM files to the user
echo "Sorted BAM files:"
ls $INPUT_DIR/*_sorted.bam

# Initialize arrays to store file names to merge and the corresponding output names.
files_to_merge_list=()
output_merged_names=()

# Prompt user to enter the list of files to be merged
while true; do
  echo "Enter the list of files to be merged (space-separated), or type 'done' to finish:"
  read -a files_to_merge
  if [[ "${files_to_merge[0]}" == "done" ]]; then
    break
  fi

  # Prompt user to enter the output name for the merged files
  echo "Enter the output name for the merged files:"
  read output_merged_name

  # Append the file list and output name to the arrays
  files_to_merge_list+=("${files_to_merge[*]}")
  output_merged_names+=("$output_merged_name")
done

# Loop through the list of files to be merged
for i in "${!files_to_merge_list[@]}"; do
  input_files_list=(${files_to_merge_list[$i]})
  output_merged_name=${output_merged_names[$i]}
  output_merged_file="${INPUT_DIR}/${output_merged_name}_pooled.bam"

  echo "Merging ${input_files_list[*]} into ${output_merged_file}"

  # Construct the input string for Picard
  input_string=""
  for input_file in "${input_files_list[@]}"; do
    input_string+="INPUT=${INPUT_DIR}/${input_file} "
  done

  # Merge the BAM files using Picard with nice and adding missing read groups
  echo "Running Picard MergeSamFiles..."
  nice java -jar $PICARD MergeSamFiles \
    $input_string \
    OUTPUT=$output_merged_file \
    SORT_ORDER=coordinate \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=LENIENT

  # Index the newly created pooled.bam file using Samtools
  echo "Running Samtools index..."
  nice $SAMTOOLS index $output_merged_file

  echo "Merging complete for ${output_merged_file}"
done
