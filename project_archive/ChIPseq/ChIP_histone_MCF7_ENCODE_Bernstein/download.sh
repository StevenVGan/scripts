#!/bin/bash

HOME="/mnt/home/digan/work"
BASE="${HOME}/seq/ChIPseq/ChIP_histone_MCF7_ENCODE_Bernstein"
INPUT_DIR="${BASE}/align/track"

# H3K27ac
wget -O ${INPUT_DIR}/ChIP_H3K27ac_MCF7.bw https://www.encodeproject.org/files/ENCFF051GGA/@@download/ENCFF051GGA.bigWig

# H3K27me3
wget -O ${INPUT_DIR}/ChIP_H3K27me3_MCF7.bw https://www.encodeproject.org/files/ENCFF952HXB/@@download/ENCFF952HXB.bigWig

# H3K36me3
wget -O ${INPUT_DIR}/ChIP_H3K36me3_MCF7.bw https://www.encodeproject.org/files/ENCFF232HHP/@@download/ENCFF232HHP.bigWig

# H3K4me1
wget -O ${INPUT_DIR}/ChIP_H3K4me1_MCF7.bw https://www.encodeproject.org/files/ENCFF407HFV/@@download/ENCFF407HFV.bigWig

# H3K4me2
wget -O ${INPUT_DIR}/ChIP_H3K4me2_MCF7.bw https://www.encodeproject.org/files/ENCFF732IHS/@@download/ENCFF732IHS.bigWig

# H3K4me3
wget -O ${INPUT_DIR}/ChIP_H3K4me3_MCF7.bw https://www.encodeproject.org/files/ENCFF312VBM/@@download/ENCFF312VBM.bigWig

# H3K9ac
wget -O ${INPUT_DIR}/ChIP_H3K9ac_MCF7.bw https://www.encodeproject.org/files/ENCFF814UFD/@@download/ENCFF814UFD.bigWig

# H3K9me2
wget -O ${INPUT_DIR}/ChIP_H3K9me2_MCF7.bw https://www.encodeproject.org/files/ENCFF134ARV/@@download/ENCFF134ARV.bigWig

# H3K9me3
wget -O ${INPUT_DIR}/ChIP_H3K9me3_MCF7.bw https://www.encodeproject.org/files/ENCFF877NQZ/@@download/ENCFF877NQZ.bigWig

# H3K36me3
wget -O ${INPUT_DIR}/ChIP_H3K36me3_MCF7.bw https://www.encodeproject.org/files/ENCFF232HHP/@@download/ENCFF232HHP.bigWig

# H4K20me1
wget -O ${INPUT_DIR}/ChIP_H4K20me1_MCF7.bw https://www.encodeproject.org/files/ENCFF220WUF/@@download/ENCFF220WUF.bigWig


# H2AFZ
wget -O ${INPUT_DIR}/ChIP_H2AFZ_MCF7.bw https://www.encodeproject.org/files/ENCFF891IUG/@@download/ENCFF891IUG.bigWig

# H3F3A
wget -O ${INPUT_DIR}/ChIP_H3F3A_MCF7.bw https://www.encodeproject.org/files/ENCFF389RCQ/@@download/ENCFF389RCQ.bigWig