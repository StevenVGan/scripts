# Pipeline

Main CUT&RUN analysis pipeline and shared utilities.

## Contents

| Path | Purpose |
|------|---------|
| **cutrun/** | Standard CUT&RUN upstream pipeline: trim → align → peak call → QC |
| **tools/** | Analysis utilities: heatmaps (`heatmap.sh`), BED liftover, peak set ops, getfasta, subsample; **go_enrichr.py** / **annotation_pie.py** (generic Enrichr + annotation pies); **chip_downstream_reference/** (frozen GSE59530 Method 2 downstream scripts) |
| **tools/prep/** | Upstream FASTQ prep: IGM FTP + ENA SRR download, Illumina lane merge, symlink into project `data/` |

See **[tools/prep/README.md](tools/prep/README.md)** for download / link workflows before **`1_trim_qc.sh`**.

## cutrun

Editable, config-driven pipeline. Copy into each project and adjust `0_config.sh`.

See [cutrun/README.md](cutrun/README.md) for step-by-step details.

## tools

Standalone scripts. Most have config blocks at the top; edit paths before running.

See [tools/README.md](tools/README.md) for usage of each tool.

## Recent tool changes

Short changelog for shared utilities lives in the parent [scripts README](../README.md#recent-pipeline-tool-updates) (`tools/prep/`).
