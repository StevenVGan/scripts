#!/bin/bash

HOME="/mnt/home/digan/work"
BASE="${HOME}/seq/ChIPseq/ChIP_SUZ12_NONO_GATAD2B_MCF7_ENCODE_Snyder"
INPUT_DIR="${BASE}/align/track"

# SUZ12
wget -O ${INPUT_DIR}/ChIP_SUZ12_MCF7.bw https://www.encodeproject.org/files/ENCFF522YXV/@@download/ENCFF522YXV.bigWig

# NONO
wget -O ${INPUT_DIR}/ChIP_NONO_MCF7.bw https://www.encodeproject.org/files/ENCFF792AJV/@@download/ENCFF792AJV.bigWig

# GATAD2B
wget -O ${INPUT_DIR}/ChIP_GATAD2B_MCF7.bw https://www.encodeproject.org/files/ENCFF380IIP/@@download/ENCFF380IIP.bigWig
