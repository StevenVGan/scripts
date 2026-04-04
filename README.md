# Analysis Scripts

Bash scripts for CUT&RUN and related sequencing analysis pipelines.

## Directory Structure

```
scripts/
├── pipeline/              # Main CUT&RUN pipeline (current standard)
│   ├── cutrun/            # Upstream: trim → align → peak call → QC
│   └── tools/             # Reusable utilities (heatmap, link_fastq, etc.)
├── docs/                  # Tutorials + GitHub Pages (rendered HTML)
├── project_archive/       # Archived project-specific scripts (by assay type)
├── legacy/                # Deprecated pipelines (MACS2-based)
└── experimental/         # Work-in-progress (peak calling tests)
```

## Quick Start

For **new CUT&RUN projects**, use the standard pipeline:

1. Copy `pipeline/cutrun/` into your project `script/` directory
2. Edit `0_config.sh` (paths, genome, sample groups)
3. Add `peakcall_groups.tsv` for peak calling
4. Run: `./run_all.sh`

## Pipeline Overview

| Folder | Purpose |
|--------|---------|
| **pipeline/cutrun** | Full CUT&RUN upstream: Trim Galore → Bowtie2 → HOMER tags → MACS3/HOMER peaks → QC + MultiQC |
| **pipeline/tools** | Standalone utilities: heatmaps, BED liftover, merge peaks, intersect peaks, link FASTQs, subsample |
| **project_archive** | One-off scripts from past projects (ChIP-seq, ATAC-seq, PRO-seq, CUT&RUN) |
| **legacy** | Old MACS2-based pipeline; kept for reference |
| **experimental** | Testing MACS3, HOMER, SEACR peak callers; peak set modes (ComplexHeatmap UpSet) |
| **docs** | Tutorials + GitHub Pages source (CUT&RUN analysis) |

## Requirements

- Conda environment with: `trim_galore`, `cutadapt`, `fastqc`, `bowtie2`, `samtools`, `deeptools`, `homer`, `bedtools`, `macs3`, `preseq`, `multiqc`, `phantompeakqualtools` (provides `run_spp.R`)
- Reference genome index (Bowtie2) and blacklist BED
- Projects expect data at `$HOME/work/seq/` (configurable in `0_config.sh`)

## Tutorials

Course tutorials in `docs/tutorials/` cover CUT&RUN analysis (AWS setup, QC, alignment, peak calling).

**Enable rendered HTML on GitHub Pages:** Settings → Pages → Source: Deploy from a branch → Branch: `main` → Folder: `/docs`. The site will be at `https://<username>.github.io/<repo>/`.

**Verify:** After pushing, visit the Pages URL and check that the CUT&RUN tutorial link loads correctly. See `docs/README.md` for local preview and troubleshooting.

## Recent pipeline tool updates

### `pipeline/tools/link_fastq.sh`

Shared helper to symlink Illumina-style or ENA-style FASTQs into a project `data/` tree from a TSV map (`prefix` TAB `newname`; optional extra columns ignored; `#` comment lines skipped).

- **Symlinks:** Uses `ln -sfn` so re-runs replace targets without leaving stale names.
- **PE (Illumina):** Resolves `${RAW_DIR}/${prefix}*_R1_001.fastq.gz` and `*_R2_001.fastq.gz` with `nullglob` so a missing match is not mistaken for a literal path (avoids broken links).
- **SE / GEO–ENA layout:** If there is no single R1 Illumina match but `${RAW_DIR}/${prefix}.fastq.gz` exists (e.g. `SRR123.fastq.gz`), that file is linked as `${newname}_R1.fastq.gz`.
- **Stale R2:** If no R2 source is found, any existing `${newname}_R2.fastq.gz` in `DEST_DIR` is removed (cleans up after PE→SE or bad earlier runs).
- **Config:** Defaults in the script can be overridden with `export RAW_DIR`, `DEST_DIR`, `MAP_FILE` so a project wrapper can `exec` this file unchanged.
- **Log:** Writes an enriched run log beside the map: `link_sample.log.tsv` (same basename as `MAP_FILE`, `.tsv` → `.log.tsv`) with resolved source paths and timestamp; the map file itself is not overwritten.

Example project pattern: `seq/<assay>/<project>/script/link_fastq.sh` sets the three variables then `exec`s this tool. ChIP-seq GSE59530 (`seq/ChIPseq/MCF7_ER_p65_ChIP_GSE59530`) uses `link_sample.tsv` and `script/link_fastq.sh`; see that project’s `RUNBOOK.txt`.

### `pipeline/tools/download_geo_fastq_ena.sh`

ENA HTTPS downloader for SRR accessions (parallel `curl`, `*.part` then rename, optional `md5sum -c`). Intended when SRA toolkit or campus FTP is awkward; **single-end / one FASTQ per SRR** in the current implementation. Override paths via `DEST_DIR`, `SRR_LIST_FILE`, `LOG_DIR`, `MD5_FILE`, `DOWNLOAD_JOBS`, or pass a list file as the first argument.

Details and copy-paste examples: [pipeline/tools/README.md](pipeline/tools/README.md).