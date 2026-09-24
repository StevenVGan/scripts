# Pipeline

Standard upstream pipelines and shared utilities. See
[../CONVENTIONS.md](../CONVENTIONS.md) for layout / file-template /
promotion-path conventions; this README is a per-pipeline contents map.

## Contents

| Path | Purpose |
|------|---------|
| **cnr/** | Standard CUT&RUN upstream pipeline: trim → align → peak call → QC; CUT&Tag and ChIP projects reuse it; see [cnr/README.md](cnr/README.md) |
| **csrna/** | Nascent/csRNA fork: post-trim MultiQC gate, `TRIM_GALORE_EXTRA` / `BT2_EXTRA`, strand bigWigs, HOMER `-sspe`; see [csrna/README.md](csrna/README.md) |
| **pro/** | PRO-seq fork: poly-A trim (1.1), strand flip, pausing index + divergent calls (4.3); see [pro/README.md](pro/README.md) |
| **atac/** | ATAC-seq pipeline (paired-end, fragment-size aware); see [atac/README.md](atac/README.md) |
| **multiome/** | Intentionally empty: single-cell / multiome work is project-local Python under `analysis/` (CONVENTIONS §12); see [multiome/README.md](multiome/README.md) |
| **templates/** | Scaffolding templates stamped by `../setup/new_project.sh` (gitignore, TSV headers, README) |
| **tools/** | Analysis utilities: heatmaps (`heatmap.sh`), BED liftover, peak set ops, getfasta, subsample; **go_enrichr.py** / **annotation_pie.py** (generic Enrichr + annotation pies); **igm_manifest.py** (IGM submission TSVs). Topic subfolders allowed for methodology families (see [../CONVENTIONS.md](../CONVENTIONS.md) §5). **chip_downstream_reference/** is now a pointer to `seq/_joint/MCF7_ER_p65_cobinding/`. |
| **tools/prep/** | Upstream FASTQ prep: IGM download via Globus (legacy FTP kept) + ENA SRR download, Illumina lane merge, symlink into project `data/` |
| **tools/viz/** | `_figure_style.py` (`apply_publication_style()`, canonical matplotlib rcParams) + `_profile_plot.py`; analyses import these or copy `_figure_style.py` verbatim into `script/` |

See **[tools/prep/README.md](tools/prep/README.md)** for download / link workflows before **`1_trim_qc.sh`**.

## Pipeline variants

All four (`cnr`, `csrna`, `pro`, `atac`) follow the same shape:
`0_config.sh` + numbered step scripts + `run_all.sh` driven by `RUN_*` env
toggles. They all source `0_config.sh` for the `bio` PATH prepend, the
`log_start` helper (with `LOG_KEEP_N=3` auto-prune), and the
`emit_references_tsv` helper (which `5_qc.sh` invokes to write
`<project>/references.tsv` per CONVENTIONS.md §4 / §9).

Per-project use: `../setup/new_project.sh <variant> <project_name>` copies
the variant's steps flat into `seq/<assay>/<project>/script/` (numbered steps
are untracked copies; `0_config.sh`, `run_all.sh`, `link_fastq.sh`, `local/`
are tracked) and bakes `GENOME`/`SE` in from its flags; `BASE` auto-derives.
See CONVENTIONS §setup and each variant's README for assay-specific knobs.
Bulk RNA-seq (STAR) has no variant here yet — CONVENTIONS §setup.

## tools

Standalone scripts. Most have config blocks at the top; edit paths before running.

See [tools/README.md](tools/README.md) for usage of each tool.
