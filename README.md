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
| **pipeline/tools** | Standalone utilities: heatmaps, BED liftover, merge peaks, link FASTQs, subsample |
| **project_archive** | One-off scripts from past projects (ChIP-seq, ATAC-seq, PRO-seq, CUT&RUN) |
| **legacy** | Old MACS2-based pipeline; kept for reference |
| **experimental** | Testing MACS3, HOMER, SEACR peak callers without control |
| **docs** | Tutorials + GitHub Pages source (CUT&RUN analysis) |

## Requirements

- Conda environment with: `trim_galore`, `cutadapt`, `fastqc`, `bowtie2`, `samtools`, `deeptools`, `homer`, `bedtools`, `macs3`, `preseq`, `multiqc`, `phantompeakqualtools` (provides `run_spp.R`)
- Reference genome index (Bowtie2) and blacklist BED
- Projects expect data at `$HOME/work/seq/` (configurable in `0_config.sh`)

## Tutorials

Course tutorials in `docs/tutorials/` cover CUT&RUN analysis (AWS setup, QC, alignment, peak calling).

**Enable rendered HTML on GitHub Pages:** Settings → Pages → Source: Deploy from a branch → Branch: `main` → Folder: `/docs`. The site will be at `https://<username>.github.io/<repo>/`.

**Verify:** After pushing, visit the Pages URL and check that the CUT&RUN tutorial link loads correctly. See `docs/README.md` for local preview and troubleshooting.