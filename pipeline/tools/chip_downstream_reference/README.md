# ChIP downstream reference (GSE59530 Method 2)

Verbatim copies of exploratory **downstream visualization** scripts from the MCF7 ER/p65 cobinding analysis. They are **not** substitutes for the originals under the project tree.

## Authoritative source

Edit and version these files in:

`seq/ChIPseq/MCF7_ER_p65_ChIP_GSE59530/analysis/p65_ER_cobinding/`

When you change behavior or fix bugs there, **refresh this folder** with the same files (manual sync). See [SOURCE.txt](SOURCE.txt) for the path manifest.

## Contents

| Subfolder | What it does |
|-----------|----------------|
| `go_enrichment/` | Enrichr GO Biological Process for `downstream/lists/{macs,homer}/genes_Method2_*.txt` |
| `annotation_composition/` | Genomic annotation composition (pies, TSVs) vs Method 2 arms |
| `signal_profiles/` | Site-class BEDs + deepTools heatmaps (rep averaging + normalization) |

## Running from this copy

Scripts use `SCRIPT_DIR` / `__file__` for sibling Python helpers, so running from this directory works if you point inputs (e.g. `--cobinding-root`, `TRACK_DIR`) at **your** data. Paths inside defaults may still mention the original machine; override with environment variables as in each script’s header comments and local `README.txt` files.

## Generic heatmaps

For track-vs-region heatmaps without the Method 2 site-class logic, use the existing tool (unchanged, self-contained):

`../heatmap.sh`

No shared `deeptools_helpers.bash` is required for this reference bundle (per standardization plan).

## Generic tools (sibling scripts)

For arbitrary inputs (not Method 2 layout), use the CLI tools in the parent `tools/` directory:

- **`../go_enrichr.py`** — gene list, annotatePeaks TSV (`Gene Name`), or BED name column → Enrichr + bar plot (`--gene-set` e.g. `GO_Biological_Process_2023`).
- **`../annotation_pie.py`** — HOMER `Annotation` column (matched case-insensitively), or another TSV column / one string per line → pie chart.

This bundle stays the frozen GSE59530 reference; those scripts are the reusable entry points for new projects.
