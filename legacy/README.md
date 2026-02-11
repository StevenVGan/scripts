# Legacy Scripts

Deprecated CUT&RUN pipeline. **Do not use for new projects.**

## What's here

- **0_link_fastq.sh**: Symlink raw FASTQs into project.
- **1_trim_bowtie2_macs2.sh**: Monolithic script combining trim, alignment, and MACS2 peak calling (FRiP).
- **2_extra_qc.sh**, **2_qc_extra.sh**: QC steps (duplicate naming).
- **downsample.sh**: Subsample FASTQs for testing.

## Why deprecated

- Uses MACS2 (migrated to MACS3 in standard pipeline).
- Single large script vs modular step scripts.
- Replaced by `pipeline/cutrun/` + `pipeline/tools/`.

## Use instead

- **Standard pipeline**: `pipeline/cutrun/` with `run_all.sh`
- **Subsampling**: `pipeline/tools/subsample_data.sh`
