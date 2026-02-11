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

## Status

Experimental. Not part of the standard pipeline. Use `pipeline/cutrun/` for production runs.
