# Tools

Reusable utilities for sequencing analysis. Each script has a config block at the top; edit paths before running.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| **heatmap.sh** | deepTools heatmaps of bigWig tracks over peaks/TSS/genes | Edit config (BASE, BED_FILES, BW_FILES, REGION), then `./heatmap.sh` |
| **lift_bed.sh** | LiftOver BED files between assemblies (e.g. hg19→hg38) | `./lift_bed.sh INPUT FROM TO OUTPUT [CHAIN]` |
| **merge_peaks.sh** | Merge multiple BED/HOMER peak files into one | `./merge_peaks.sh OUTPUT.bed INPUT1 [INPUT2 ...]` |
| **intersect_peaks.sh** | Intersect peaks (regions in ALL inputs) + UpSet plot | `./intersect_peaks.sh [--slop BP] [--names "N1,N2,..."] OUTPUT.bed INPUT1 INPUT2 [...]` |
| **peak_ops.sh** | Peak set ops: intersect, distinct, or union + UpSet plot | `./peak_ops.sh --mode MODE [--slop BP] [--names "N1,N2,..."] OUTPUT.bed INPUT1 INPUT2 [...]` |
| **getfasta.sh** | Extract sequences from reference genome for BED regions | `./getfasta.sh [--slop BP] [--no-center] [--tab] OUTPUT.txt INPUT1.bed [INPUT2.bed ...]` |
| **link_fastq.sh** | Symlink/copy raw FASTQs into project `data/` with sample mapping | Edit config (RAW_DIR, DEST_DIR, MAP_FILE), then `./link_fastq.sh` |
| **subsample_data.sh** | Subsample FASTQs (e.g. 1M reads) for testing | Uses `0_config.sh` from `../cutrun/` or set `CONFIG_FILE=/path/to/project/0_config.sh`; needs `seqtk` |

## Examples

**Heatmap**
```bash
# Edit BASE, BED_FILES, BW_FILES, REGION (Peaks|TSS|Genes), OUTPUT_NAME
./heatmap.sh
```

**Lift BED**
```bash
./lift_bed.sh ./MCF7_Amir_hg19 hg19 hg38 ./MCF7_Amir_hg38
./lift_bed.sh my.bed hg19 hg38 ./out
```

**Merge peaks**
```bash
./merge_peaks.sh merged.bed rep1_peaks.bed rep2.annotatePeaks.txt
```

**Intersect peaks** (regions in ALL inputs; generates UpSet PDF)
```bash
./intersect_peaks.sh out.bed rep1.annotatePeaks.txt rep2.annotatePeaks.txt
./intersect_peaks.sh --slop 250 --names "R1,R2,R3,R4" peaks/ERa_intersect.bed peaks/*.annotatePeaks.txt
```

**Peak set operations** (intersect / distinct / union; generates UpSet PDF)
```bash
# intersect: regions in ALL inputs
./peak_ops.sh --mode intersect --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/intersect.bed peaks/*.annotatePeaks.txt

# distinct: regions in exactly the specified sets (partition: in A&B, NOT in C or D)
./peak_ops.sh --mode distinct --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/distinct.bed peaks/*.annotatePeaks.txt

# union: regions in ANY input
./peak_ops.sh --mode union --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/union.bed peaks/*.annotatePeaks.txt
```
Requires: R + UpSetR (`install.packages('UpSetR', repos='https://cloud.r-project.org')`)

**Get sequences from BED regions**
```bash
# Default: center regions, extend ±500 bp (1001 bp total), output sequences only
./getfasta.sh seqs.txt peaks.bed

# Use original coordinates, no slop
./getfasta.sh --no-center --slop 0 seqs.txt peaks.bed

# Tab format (name, seq) for downstream use
./getfasta.sh --tab seqs.txt rep1.bed rep2.bed
```

**Link FASTQs**
```bash
# Edit RAW_DIR, DEST_DIR, MAP_FILE (format: prefix TAB newname TAB ...)
./link_fastq.sh
```

## Notes

- **link_fastq** and **subsample_data**: project-specific; edit config paths before running.
- **subsample_data**: Uses `../cutrun/0_config.sh` when run from `pipeline/tools/`, or set `CONFIG_FILE=/path/to/project/script/0_config.sh` for your project.
