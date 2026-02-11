#!/bin/bash

# Description: This script takes all the _pooled_filtered.bed files and asks the user for each file whether they want to intersect it with other files in searching for high confidence set.
# Input: Filtered BED file with names ending *_filtered.bed
# Last modified time: 2024-04-04

# Set directories and parameters
HOME="/mnt/home/digan/work"
BASE="${HOME}/seq/CUTRUN/240202_240227_CUTRUN_NCAPH2_HA_HCT116_WT_AID_Steven"
INPUT_DIR="${BASE}/peaks"
BEDTOOLS="/opt/apps/bio/bedtools-2.27.1/bin"
OVERLAPRATE=0.3

# Set the log file name with the output directory path
LOGFILE="$INPUT_DIR/4_hiConfi_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

# Find and Print all the _pooled_filtered.bed files and _indiv_filtered.bed files in the input directory
pooled_bed_files=($(find $INPUT_DIR -type f -name "*_pooled_filtered.bed"))
overlap_files=($(find $INPUT_DIR -type f -name "*_indiv_filtered.bed"))
echo "The following pooled bed files were found:"
ls $INPUT_DIR/*_pooled_filtered.bed
echo "The following individual bed files were found:"
ls $INPUT_DIR/*_indiv_filtered.bed

# For each _pooled_filtered.bed file, ask the user whether they want to intersect (overlap) it with other files
for pooled_bed in "${pooled_bed_files[@]}"; do
  pooled_bed="$(basename ${pooled_bed})"
  echo "Do you want to intersect (overlap) ${pooled_bed} with other files? (yes/no)"
  read user_decision

  if [ "$user_decision" == "yes" ]; then
    # Ask the user to enter the files they want to intersect with, space-separated
    echo "Enter the file names you want to intersect (overlap) with ${pooled_bed}, space-separated:"
    read -a selected_overlap_files

    # Perform the intersection for each specified file
    for overlap_file in "${selected_overlap_files[@]}"; do
      output_name=$(basename "${pooled_bed%_pooled_filtered.bed}")
      output_file="${INPUT_DIR}/${output_name}_hiconfi.bed"
      nice $BEDTOOLS/intersectBed -a ${INPUT_DIR}/$pooled_bed -b ${INPUT_DIR}/$overlap_file -wa -f $OVERLAPRATE -r > $output_file
      echo "Intersected ${pooled_bed} with ${overlap_file}, output: ${output_file}"
    done
  else
    echo "No intersection will be performed for ${pooled_bed}."
  fi
done

echo "Condifent peaks searching complete!"
