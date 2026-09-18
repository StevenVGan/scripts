# Analysis Scripts

Bash + Python scripts for CUT&RUN / CUT&Tag / ChIP-seq, ATAC, csRNA and
PRO-seq upstream pipelines, plus shared tools, project scaffolding, conda-env
tracking and maintenance used across `~/work/seq/` projects.

## Conventions (read this first)

For all "where does X go?" / "how do I structure Y?" questions, see
**[CONVENTIONS.md](CONVENTIONS.md)** — the single source of truth for
`~/work/` layout, gitignore rules, file templates, promotion path,
joint-analysis conventions, env/ref tracking, and log policy. Both this
README and per-project READMEs assume CONVENTIONS.md.

## Directory Structure

```
scripts/
├── CONVENTIONS.md         # source of truth for ~/work/ layout (read first)
├── env/                   # conda env tracking: bio rna sc meth primer pwm2 pwm2-pb vienna + dated lock/ — see env/README.md
├── setup/                 # new_project.sh scaffolder, gitall.sh, check_env_drift.sh, clone_all.sh, sync_brain.sh, bootstrap.md
├── maintenance/           # archive_inactive.sh (BAM→CRAM for inactive projects), dedup/dupescan.py
├── pipeline/              # All standard pipelines + shared tools
│   ├── cnr/               # CUT&RUN (also CUT&Tag / ChIP via --dest): trim → align → peak call → QC
│   ├── csrna/             # csRNA fork: post-trim MultiQC gate, strand bigWigs
│   ├── pro/               # PRO-seq fork: strand flip, pausing index
│   ├── atac/              # ATAC-seq pipeline
│   ├── multiome/          # intentionally empty scaffold (single-cell work is project-local; CONVENTIONS §12)
│   ├── templates/         # canonical scaffolding templates (gitignore, TSVs, README)
│   └── tools/             # Reusable utilities; prep/ + viz/ subfolders; topic subfolders OK (see CONVENTIONS.md §5)
├── docs/                  # Tutorials + GitHub Pages (rendered HTML)
├── project_archive/       # Archived project-specific scripts (by assay type)
├── legacy/                # Deprecated pipelines (MACS2-based)
└── experimental/          # Work-in-progress (peak calling tests)
```

## Quick Start

For a **new sequencing project**, see CONVENTIONS.md `§setup`. The short
version:

1. `scripts/setup/new_project.sh <cnr|csrna|pro|atac> <project_name>` (`--dest cnt|chip`, `--se`, `--genome`)
   — stamps the flat `script/` (numbered steps as untracked copies), `.gitignore`, template TSVs, README
2. Author `link_sample.tsv`, `samples.tsv` and `peakcall_groups.tsv`; set `RAW_DIR` in `script/link_fastq.sh`
3. `cd script && ./link_fastq.sh && ./run_all.sh` — `references.tsv` is auto-emitted by `5_qc.sh`
4. Bulk RNA-seq has no pipeline dir yet (CONVENTIONS §setup); single-cell / re-analysis projects are built by hand (§12, §setup)

## Pipeline Overview

| Folder | Purpose |
|--------|---------|
| **pipeline/cnr** | Full CUT&RUN upstream: Trim Galore → Bowtie2 → HOMER tags → MACS3/HOMER peaks → QC + MultiQC. CUT&Tag and ChIP projects reuse it (`new_project.sh cnr <name> --dest cnt`) |
| **pipeline/csrna** | csRNA fork of cnr: post-trim MultiQC gate, strand bigWigs, HOMER `-sspe` |
| **pipeline/pro** | PRO-seq fork: poly-A trim (1.1), strand flip, pausing index (4.3), divergent calls |
| **pipeline/atac** | ATAC-seq pipeline (paired-end, fragment-size aware) |
| **pipeline/templates** | Scaffolding templates consumed by `setup/new_project.sh` (gitignore, `samples` / `link_sample` / `peakcall_groups` / `tss_groups` TSV headers, project README; `README.analysis.md` is copied by hand per CONVENTIONS §11) |
| **pipeline/tools** | Standalone utilities: heatmaps, BED liftover, peak set ops, getfasta, **go_enrichr.py** / **annotation_pie.py**, **igm_manifest.py** (IGM submission TSVs); **prep/** = IGM/ENA FASTQ download + link/merge; **viz/** = shared matplotlib style + profile plotter; topic subfolders for methodology families per CONVENTIONS.md §5 |
| **setup/** | `new_project.sh` (project scaffolder), `gitall.sh` (dirty/ahead report for repos up to two levels below `~/work`; per-project repos under `seq/<assay>/<project>/` are not scanned), `check_env_drift.sh` (env newer than its yml), `clone_all.sh` + `bootstrap.md` (second-node setup), `sync_brain.sh` |
| **maintenance/** | `archive_inactive.sh` (lossless BAM→CRAM + cleanup for `projects.tsv` rows marked `cleanup`), `dedup/dupescan.py` (cross-volume duplicate survey) |
| **env/** | One yml per tracked env (`bio`, `rna`, `sc`, `meth`, `primer`, `pwm2`, `pwm2-pb`, `vienna`) + `lock/<env>.<date>.yml` snapshots. Spec vs export rules in `env/README.md` |
| **project_archive** | One-off scripts from past projects (ChIP-seq, ATAC-seq, PRO-seq, CUT&RUN) |
| **legacy** | Old MACS2-based pipeline; kept for reference |
| **experimental** | Testing MACS3, HOMER, SEACR peak callers; peak set modes (ComplexHeatmap UpSet) |
| **docs** | Tutorials + GitHub Pages source (CUT&RUN analysis) |

## Requirements

- The **`bio` conda env** — the env the bulk pipelines and tools use
  (`trim_galore`, `cutadapt`, `fastqc`, `bowtie2`, `samtools`,
  `deeptools`, `homer`, `bedtools`, `macs3`, `preseq`, `multiqc`,
  `phantompeakqualtools` providing `run_spp.R`, plus R packages).
  Tracked in `env/bio.yml`; see `env/README.md` for rebuild and
  versioning instructions (and for the other tracked envs: `rna` for
  STAR RNA-seq, `sc` for single-cell, `meth`, `primer`, `pwm2*`, `vienna`).
  Per-project `references.tsv` (auto-emitted by `5_qc.sh`) records which
  env lock was active at run time. `0_config.sh` prepends the env's
  `bin/` to `PATH`; nothing `conda activate`s.
- Reference genome index (Bowtie2) and blacklist BED — paths in
  `0_config.sh`. Per-project `references.tsv` records which were used.
- Projects expect data at `$HOME/work/seq/` (configurable in
  `0_config.sh`).

## Tutorials

Course tutorials in `docs/tutorials/` cover CUT&RUN analysis (AWS setup, QC, alignment, peak calling).

**Enable rendered HTML on GitHub Pages:** Settings → Pages → Source: Deploy from a branch → Branch: `main` → Folder: `/docs`. The site will be at `https://<username>.github.io/<repo>/`.

**Verify:** After pushing, visit the Pages URL and check that the CUT&RUN tutorial link loads correctly. See `docs/README.md` for local preview and troubleshooting.

## Recent pipeline tool updates

### `pipeline/tools/prep/` (FASTQ prep)

Upstream helpers live under **`pipeline/tools/prep/`**: **`download_fastq.sh`** (IGM FTP), **`download_geo_fastq_ena.sh`** (ENA SRR / SE), **`link_fastq.sh`** (Illumina `*_R1_001` / SE maps), **`merge_lanes_inplace.sh`**, **`link_merged_fastqs.sh`**. See [pipeline/tools/prep/README.md](pipeline/tools/prep/README.md).

### `pipeline/tools/prep/link_fastq.sh`

Shared helper to symlink Illumina-style or ENA-style FASTQs into a project `data/` tree from a TSV map (`prefix` TAB `newname`; optional extra columns ignored; `#` comment lines skipped).

- **Symlinks:** Uses `ln -sfn` so re-runs replace targets without leaving stale names.
- **PE (Illumina):** Resolves `${RAW_DIR}/${prefix}*_R1_001.fastq.gz` and `*_R2_001.fastq.gz` with `nullglob` so a missing match is not mistaken for a literal path (avoids broken links).
- **SE / GEO–ENA layout:** If there is no single R1 Illumina match but `${RAW_DIR}/${prefix}.fastq.gz` exists (e.g. `SRR123.fastq.gz`), that file is linked as `${newname}_R1.fastq.gz`.
- **Stale R2:** If no R2 source is found, any existing `${newname}_R2.fastq.gz` in `DEST_DIR` is removed (cleans up after PE→SE or bad earlier runs).
- **Config:** Defaults in the script can be overridden with `export RAW_DIR`, `DEST_DIR`, `MAP_FILE` so a project wrapper can `exec` this file unchanged.
- **Log:** Writes an enriched run log beside the map: `link_sample.log.tsv` (same basename as `MAP_FILE`, `.tsv` → `.log.tsv`) with resolved source paths and timestamp; the map file itself is not overwritten.

Example project pattern: `seq/<assay>/<project>/script/link_fastq.sh` sets the three variables then `exec`s **`prep/link_fastq.sh`**. ChIP-seq GSE59530 (`seq/chip/ChIP_GSE59530_ER_p65_MCF7`) uses `link_sample.tsv` and `script/link_fastq.sh`; see that project’s `RUNBOOK.txt`.

### `pipeline/tools/prep/download_geo_fastq_ena.sh`

ENA HTTPS downloader for SRR accessions (parallel `curl`, `*.part` then rename, optional `md5sum -c`). Intended when SRA toolkit or campus FTP is awkward; **single-end / one FASTQ per SRR** in the current implementation. Override paths via `DEST_DIR`, `SRR_LIST_FILE`, `LOG_DIR`, `MD5_FILE`, `DOWNLOAD_JOBS`, or pass a list file as the first argument.

Details and copy-paste examples: [pipeline/tools/README.md](pipeline/tools/README.md).

### `pipeline/tools/intersect_peaks.sh`

Thin wrapper around [`peak_ops.sh --mode intersect`](pipeline/tools/peak_ops.sh): same BED output, UpSet, and optional Venn (`--viz`). Prefer calling `peak_ops.sh` directly when you need `--mode distinct` or `union`.

### `pipeline/tools/go_enrichr.py` and `annotation_pie.py`

Generic **Enrichr** enrichment (gene list, HOMER annotate table, or BED name column) plus a horizontal bar plot of top terms, and **annotation composition pies** from HOMER-style `Annotation` strings (or an arbitrary TSV column). See [pipeline/tools/README.md](pipeline/tools/README.md) for examples.

### `pipeline/tools/chip_downstream_reference/` — moved

The frozen ChIP downstream scripts that previously lived here have been
carved out into the dedicated joint analysis repo at
`seq/_joint/MCF7_ER_p65_cobinding/` (single source of truth — no more
edit-here-copy-there sync). This dir keeps a one-page README pointer for
discoverability. See [pipeline/tools/chip_downstream_reference/README.md](pipeline/tools/chip_downstream_reference/README.md).