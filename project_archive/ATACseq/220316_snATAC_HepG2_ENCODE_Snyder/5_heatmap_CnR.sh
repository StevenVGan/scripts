#!/bin/bash

# Description: This script takes bigWig files and compute matrix for heatmap and histogram drawing.
# Input: bigWig files with name endning in *.bw
# Last modified time: 2023-06-09

# Set directories
HOME="/mnt/home/digan"
BASE="${HOME}/seq/ATACseq/220316_snATAC_HepG2_ENCODE_Snyder"
INPUT_DIR="${BASE}/align/track"
INPUT_DIR2="${BASE}/peaks"
OUTPUT_DIR="${BASE}/matrix"
DEEPTOOLS="/opt/apps/bio/deeptools-3.5.1/bin"

# Set parameters
PEAK_UP=1000
PEAK_DN=1000
TSS_UP=1000
TSS_DN=1000
GENES_UP=3000
GENES_DN=3000
BODY_LENGTH=7500
BIN_SIZE=10
# GENOME_BED="${HOME}/ref/annot/hg19/hg19_annot_genome.bed"
GENOME_BED="${HOME}/ref/annot/hg38/hg38_annot_genome.bed"
CPU=8
COLOR="RdBu"
SORT="descend"
ZMIN=30

# Creat the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/heatmap
mkdir -p $OUTPUT_DIR/sort

# Prompt user to enter the list of files to analyse
echo "bigWig files found in the input ditectory:"
ls $INPUT_DIR/*.bw
echo "Enter the list of files to be analysed (space-separated):"
read -a bw_names

# Prompt user to enter the label of files for figure
echo "Enter the labels of these files (space-separated):"
read -a bw_labels

# Prompt user to enter the region to plot
echo "Plot what region? (Peaks/TSS/Genes/skip):"
read region

REGIONSLABEL=""
# Compute matrix for each case
if [ "$region" == "Peaks" ]; then

#  SORT="no"

  # Prompt user to enter the bed files used to compute matrix
  echo "bed files found in the input directory 2:"
  ls $INPUT_DIR2/*.bed $OUTPUT_DIR/*.bed
  echo "Enter the bed files used to compute matrix:"
  read -a bed_names
  echo "Enter the bed file labels:"
  read -a bed_labels


  # Prompt user to enter the output name
  echo "Enthr the file name of output:"
  read output_name
  output_file="${OUTPUT_DIR}/${output_name}_mtx.gz"
  output_sort="${OUTPUT_DIR}/sort/${output_name}_sorted.bed"

  bw_files=""
  for bw_name in "${bw_names[@]}"; do
    bw_files+="${INPUT_DIR}/${bw_name} "
  done

  bw_labelss=""
  for bw_label in "${bw_labels[@]}"; do
    bw_labelss+="${bw_label} "
  done

  bed_files=""
  for bed_name in "${bed_names[@]}"; do
    bed_files+="${INPUT_DIR2}/${bed_name} "
  done

  bed_labelss=""
  for bed_label in "${bed_labels[@]}"; do
    bed_labelss+="${bed_label} "
  done

  REGIONSLABEL="--regionsLabel ${bed_labelss}"

  echo "Running compute matrix on designated bed fils ..."

  nice $DEEPTOOLS/computeMatrix reference-point --referencePoint "center" -b $PEAK_UP -a $PEAK_DN \
    -S $bw_files -R $bed_files -o $output_file --outFileSortedRegions $output_sort \
    --missingDataAsZero -p $CPU --samplesLabel $bw_labelss -bs $BIN_SIZE

elif [ "$region" == "TSS" ]; then

  # Prompt user to enter the output name
  echo "Enthr the file name of output:"
  read output_name
  output_file="${OUTPUT_DIR}/${output_name}_TSS_mtx.gz"
  output_sort="${OUTPUT_DIR}/sort/${output_name}_TSS_sorted.bed"

  bw_files=""
  for bw_name in "${bw_names[@]}"; do
    bw_files+="${INPUT_DIR}/${bw_name} "
  done

  bw_labelss=""
  for bw_label in "${bw_labels[@]}"; do
    bw_labelss+="${bw_label} "
  done

  echo "Running compute matrix on TSS regions ..."

  nice $DEEPTOOLS/computeMatrix reference-point --referencePoint "TSS" -b $TSS_UP -a $TSS_DN \
    -S $bw_files -R $GENOME_BED -o $output_file --outFileSortedRegions $output_sort \
    --missingDataAsZero -p $CPU --samplesLabel $bw_labelss -bs $BIN_SIZE

elif [ "$region" == "Genes" ]; then

  # Prompt user to enter the output name
  echo "Enthr the file name of output:"
  read output_name
  output_file="${OUTPUT_DIR}/${output_name}_Genes_mtx.gz"
  output_sort="${OUTPUT_DIR}/sort/${output_name}_Genes_sorted.bed"

  bw_files=""
  for bw_name in "${bw_names[@]}"; do
    bw_files+="${INPUT_DIR}/${bw_name} "
  done

  bw_labelss=""
  for bw_label in "${bw_labels[@]}"; do
    bw_labelss+="${bw_label} "
  done

  echo "Running compute matrix on Genes regions ..."

  nice $DEEPTOOLS/computeMatrix scale-regions -m $BODY_LENGTH -b $GENES_UP -a $GENES_DN \
    --startLabel "TSS" --endLabel "TES" \
    -S $bw_files -R $GENOME_BED -o $output_file --outFileSortedRegions $output_sort \
    --missingDataAsZero -p $CPU --samplesLabel $bw_labelss -bs $BIN_SIZE

else
  echo "Error! Skipping computeMatrix"
fi

echo "computeMatrix complete"

# Prompt user to enter the list of files to plot
echo "_mtx.gz files found in the input ditectory:"
ls $OUTPUT_DIR/*.gz
echo "Enter the list of files to be analysed (space-separated):"
read -a mtx_names

for mtx_name in "${mtx_names[@]}"; do
  mtx_file="${OUTPUT_DIR}/$mtx_name"
  plot_out="${OUTPUT_DIR}/heatmap/${mtx_name}.png"

  echo "Plotting heatmap for ${mtx_name} ..."
  nice $DEEPTOOLS/plotHeatmap -m $mtx_file -out $plot_out --colorMap $COLOR \
    --zMin -$ZMIN --zMax $ZMIN --missingDataColor 1 $REGIONSLABEL --sortRegions $SORT

done

echo "All analyses complete!"
