# Experimental / Peak Calling Tests

Work-in-progress scripts for testing peak callers and workflows.

## Purpose

Compare and test peak calling without control samples:

- **MACS3** (woCtrl, noctrl)
- **HOMER** (woCtrl)
- **SEACR** (woCtrl)

## Scripts

| Script | Description |
|--------|-------------|
| `0_link_fastq.sh` | Link raw FASTQs into project (edit config) |
| `1_trim_bowtie2_macs3.sh` | Trim Galore + Bowtie2 + MACS3 |
| `2_qc_extra.sh` | Extra QC |
| `3_peakcall_*.sh` | Peak calling variants (MACS3, HOMER, SEACR, no control) |
| `peaks_modes.sh` | Peak set ops (distinct/intersect/union) + ComplexHeatmap UpSet; experimental |

**peaks_modes** (experimental)
```bash
./peaks_modes.sh --mode intersect --slop 250 --names "R1,R2,R3,R4" peaks/out peaks/*.bed
./peaks_modes.sh --mode union peaks/out peaks/*.bed
```
Requires: R + ComplexHeatmap + GenomicRanges (`BiocManager::install(c('ComplexHeatmap','GenomicRanges'))`)

## Status

Experimental. Not part of the standard pipeline. Use `pipeline/cutrun/` for production runs.
