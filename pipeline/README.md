# Pipeline

Standard upstream pipelines and shared utilities. See
[../CONVENTIONS.md](../CONVENTIONS.md) for layout / file-template /
promotion-path conventions; this README is a per-pipeline contents map.

## Contents

| Path | Purpose |
|------|---------|
| **cutrun/** | Standard CUT&RUN upstream pipeline: trim → align → peak call → QC |
| **csRNA/** | Nascent/csRNA fork: post-trim MultiQC gate, `TRIM_GALORE_EXTRA` / `BT2_EXTRA`, strand bigWigs, HOMER `-sspe`; see [csRNA/README.md](csRNA/README.md) |
| **proseq/** | PRO-seq fork: poly-A trim (1.1), strand flip, pausing index + divergent calls (4.3); see [proseq/README.md](proseq/README.md) |
| **atac/** | ATAC-seq pipeline (paired-end, fragment-size aware); see [atac/README.md](atac/README.md) |
| **tools/** | Analysis utilities: heatmaps (`heatmap.sh`), BED liftover, peak set ops, getfasta, subsample; **go_enrichr.py** / **annotation_pie.py** (generic Enrichr + annotation pies). Topic subfolders allowed for methodology families (see [../CONVENTIONS.md](../CONVENTIONS.md) §5). **chip_downstream_reference/** is now a pointer to `seq/_joint/MCF7_ER_p65_cobinding/`. |
| **tools/prep/** | Upstream FASTQ prep: IGM FTP + ENA SRR download, Illumina lane merge, symlink into project `data/` |

See **[tools/prep/README.md](tools/prep/README.md)** for download / link workflows before **`1_trim_qc.sh`**.

## Pipeline variants

All four (`cutrun`, `csRNA`, `proseq`, `atac`) follow the same shape:
`0_config.sh` + numbered step scripts + `run_all.sh` driven by `RUN_*` env
toggles. They all source `0_config.sh` for the `bio` PATH prepend, the
`log_start` helper (with `LOG_KEEP_N=3` auto-prune), and the
`emit_references_tsv` helper (which `5_qc.sh` invokes to write
`<project>/references.tsv` per CONVENTIONS.md §4 / §9).

Per-project use: copy the variant's contents into
`seq/<assay>/<project>/script/pipeline/`, edit the project's
`0_config.sh`. See each variant's README for assay-specific knobs.

## tools

Standalone scripts. Most have config blocks at the top; edit paths before running.

See [tools/README.md](tools/README.md) for usage of each tool.
