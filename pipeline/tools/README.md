# Tools

Reusable utilities for sequencing analysis. Each script has a config block at the top; edit paths before running.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| **heatmap.sh** | deepTools heatmaps of bigWig tracks over peaks/TSS/genes | Edit config (BASE, BED_FILES, BW_FILES, REGION), then `./heatmap.sh` |
| **lift_bed.sh** | LiftOver BED files between assemblies (e.g. hg19→hg38) | `./lift_bed.sh INPUT FROM TO OUTPUT [CHAIN]` |
| **merge_peaks.sh** | Merge multiple BED/HOMER peak files into one | `./merge_peaks.sh OUTPUT.bed INPUT1 [INPUT2 ...]` |
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

**Link FASTQs**
```bash
# Edit RAW_DIR, DEST_DIR, MAP_FILE (format: prefix TAB newname TAB ...)
./link_fastq.sh
```

## Notes

- **link_fastq** and **subsample_data**: project-specific; edit config paths before running.
- **subsample_data**: Uses `../cutrun/0_config.sh` when run from `pipeline/tools/`, or set `CONFIG_FILE=/path/to/project/script/0_config.sh` for your project.
