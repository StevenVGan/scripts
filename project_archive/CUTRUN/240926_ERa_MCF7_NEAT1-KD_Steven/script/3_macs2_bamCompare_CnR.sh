#!/bin/bash

# Description: This script performs MACS2 analysis on user-defined IP data groups and their corresponding control groups, and filter the peaks. Then it will perform bamCompare upon user's decision.
# Input: Pooled or sorted BAM file with names ending *_pooled.bam or *sorted_bam
# Last modified time: 2024-10-22

# Set directories
HOME="/mnt/home/digan/work"
BASE="${HOME}/seq/CUTRUN/240926_ERa_MCF7_NEAT1-KD_Steven"
INPUT_DIR="${BASE}/align"
OUTPUT_DIR="${BASE}/peaks"
OUTPUT_DIR2="${BASE}/align/track"
MACS2="/opt/apps/bio/macs2-2.1.2/bin/macs2"
BEDTOOLS="/opt/apps/bio/bedtools-2.27.1/bin"
DEEPTOOLS="/opt/apps/bio/deeptools-3.5.1/bin"

# Set parameterers for macs2
GENOMESIZE="hs"
FORMAT="BAMPE"
FDR=0.05
BLACKLIST="${HOME}/ref/blacklist/h19/wgEncodeHg19ConsensusSignalArtifactRegions.bed"

# Set parameters for bamCompare
SAMPLE_LENGTH=1000
SAMPLE_NUMBERS=100000
BIN_SIZE=10
SMOOTH_LENGTH=30
GENOME_SIZE=2864785220
IGNORE="chrM"
CPU=16

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR2

# Set the log file name with the output directory path
LOGFILE="$OUTPUT_DIR/3_macs2_bamCompare_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

# Ask user whether to perform bamCompare or not
echo "Perform bamCompare? (yes/no)"
read user_decision

# Print all pooled.bam and sorted.bam files to the user
echo "These pooled and sorted BAM files were found in the input directory:"
ls $INPUT_DIR/*_sorted.bam $INPUT_DIR/*_pooled.bam

# Initialize arrays to store IP data groups and their corresponding control groups.
ip_data_groups=()
control_groups=()

# Prompt user to enter the IP data groups and their corresponding control groups
while true; do
  echo "Enter the IP data group (space-separated) or type 'done' to finish:"
  read -a ip_data_group
  if [[ "${ip_data_group[0]}" == "done" ]]; then
    break
  fi

  # Prompt user to enter the control group for the IP data group
  echo "Enter the control group for the IP data group:"
  read control_group

  # Append the IP data group and control group to the arrays
  ip_data_groups+=("${ip_data_group[*]}")
  control_groups+=("$control_group")
done

# Loop through the IP data groups and their corresponding control groups
for i in "${!ip_data_groups[@]}"; do
  ip_data_group_list=(${ip_data_groups[$i]})
  control_group=${control_groups[$i]}

  # Perform MACS2 analysis on each IP data group with its corresponding control group
  for ip_data_group in "${ip_data_group_list[@]}"; do
    ip_file="${INPUT_DIR}/${ip_data_group}"
    control_file="${INPUT_DIR}/${control_group}"

    output_name="${ip_data_group/_sorted.bam/_indiv}"
    output_name="${output_name/_pooled.bam/_pooled}"

    echo "Running MACS2 analysis on ${ip_data_group} with control ${control_group}..."

    # Run MACS2
    nice $MACS2 callpeak -t $ip_file -c $control_file -g $GENOMESIZE -q $FDR --format $FORMAT \
      --outdir $OUTPUT_DIR --name $output_name --call-summits

    echo "Filtering on ${OUTPUT_DIR}/${output_name}..."

    # Filter narrowPeaks using bedtools intersect
    narrowPeak_file="${OUTPUT_DIR}/${output_name}_peaks.narrowPeak"
    filtered_bed="${OUTPUT_DIR}/${output_name}_filtered.bed"

    # Run Bedtools Intersect to filter out blacklist
    nice $BEDTOOLS/intersectBed -a $narrowPeak_file -b $BLACKLIST -v > $filtered_bed

    echo "MACS2 analysis and filtering complete for ${ip_data_file} with control ${control_group}"
  done
done

echo "Peak calling complete!"

if [ "$user_decision" == "yes" ]; then
# Loop through the IP data groups and their corresponding control groups
  for i in "${!ip_data_groups[@]}"; do
    ip_data_group_list=(${ip_data_groups[$i]})
    control_group=${control_groups[$i]}

    # Perform bamCompare on each IP data group with its corresponding control group
    for ip_data_group in "${ip_data_group_list[@]}"; do
      ip_file="${INPUT_DIR}/${ip_data_group}"
      control_file="${INPUT_DIR}/${control_group}"

      bamCom_name="${ip_data_group/_sorted.bam/_indiv_cf.bw}"
      bamCom_name="${bamCom_name/_pooled.bam/_pooled_cf.bw}"
      bamCom_file="${OUTPUT_DIR2}/$bamCom_name"

      echo "Running bamCompare for ${ip_data_group} with control ${control_group}"

      nice $DEEPTOOLS/bamCompare -b1 $ip_file -b2 $control_file -o $bamCom_file -ignore $IGNORE \
        --scaleFactorsMethod "SES" -l $SAMPLE_LENGTH -n $SAMPLE_NUMBERS -bs $BIN_SIZE -p $CPU \
        --effectiveGenomeSize $GENOME_SIZE --smoothLength $SMOOTH_LENGTH --centerReads

    done
  done
  echo "bamCompare for all samples complete!"
fi
