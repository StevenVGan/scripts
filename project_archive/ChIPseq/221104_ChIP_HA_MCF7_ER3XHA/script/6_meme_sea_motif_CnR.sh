#!/bin/bash

# Description: This script asks users to select .bed files for motif analysis, extract their genomic sequences and perform meme and seacr motif analyses on them.
# Input: BED files with the name eniding *.bed
# Last modified time: 2023-04-18

# Set directories
HOME="/mnt/home/digan"
BASE="${HOME}/seq/CUTRUN/221220_230207_CUTRUN_HA_HepG2_ER3XHA_Steven"
INPUT_DIR="${BASE}/peaks"
OUTPUT_DIR="${BASE}/motif"
BEDTOOLS="/opt/apps/bio/bedtools-2.27.1/bin"
MEME_CHIP="/opt/apps/bio/meme-5.3.3/bin/meme-chip"
SEACR="/opt/apps/bio/seacr/v1.3/SEACR_1.3.sh"

# Set parameters
GENOME=hg19
GENOME_FA="/opt/archive/ref/genome/$GENOME/$GENOME.fa"
NMOTIFS=10
MAXW=30
MINW=6
MOD="zoops"
MOTIFDB="${HOME}/ref/motifs_meme/human/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
CPU=16

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Set the log file name with the output directory path
LOGFILE="$OUTPUT_DIR/6_meme_sea_motif_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

# Print all .bed files in the input directory
echo "The following bed files were found:"
ls $INPUT_DIR/*.bed

# Ask user to enter the files they want to perform motif searching on, space-seperated
echo "Enter the file names you want to perform motif searching on, space-seperated"
read -a bed_files

# Ask user whether or not to leave the fasta file undeleted
echo "Delete the fasta file (yes/no)"
read user_decision

if [ "$user_decision" == "yes" ]; then
  result=true
elif [ "$user_decision" == "no" ]; then
  result=false
else
  echo "Invalid"
  exit 1
fi

# Use for loop to process all the bed files
for bed in "${bed_files[@]}"; do
  out_name="${bed/.bed/}"
  fa_out="${OUTPUT_DIR}/${out_name}.fa.out"
  meme_out="${OUTPUT_DIR}/${out_name}_meme"
  bed="${INPUT_DIR}/$bed"

  # Run getfasta
  echo "Running getfasta for sample ${out_name} ..."
  nice $BEDTOOLS/bedtools getfasta -fi $GENOME_FA -bed $bed -fo $fa_out

  # Run meme-chip
  echo "Running meme-chip for sample ${out_name} ..."
  nice $MEME_CHIP $fa_out -meme-nmotifs $NMOTIFS -maxw $MAXW -minw $MINW -meme-mod $MOD \
    -centrimo-flip -db $MOTIFDB -o $meme_out -meme-p $CPU -minw 8

  # Delete fasta files
  if [ $result ]; then
    rm $fa_out
  fi
done

echo "All analyses complete!"
