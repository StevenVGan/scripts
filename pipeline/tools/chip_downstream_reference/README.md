# ChIP downstream reference — moved

This directory previously held frozen copies of Method 2 downstream
visualization scripts (annotation_composition / go_enrichment /
signal_profiles) from the MCF7 ER/p65 cobinding analysis.

## New location

Everything has moved to the dedicated joint analysis repo:

`~/work/seq/_joint/MCF7_ER_p65_cobinding/`

That repo is the **single source of truth** for the cobinding analysis
scripts and outputs. No more "edit-here, copy-there" sync — just edit
in the joint repo.

## When

Moved 2026-04-28 during Phase B of the ~/work restructure. See
`~/work/RESTRUCTURE_PLAN.md` §5.1 (Joint #2) and the joint repo's
`README.md` for full context.

## Generic downstream tools

Tools that aren't specific to GSE59530 still live alongside this folder:

- `../heatmap.sh` — generic deepTools bigWig heatmaps
- `../go_enrichr.py` — generic Enrichr GO enrichment from gene lists
- `../annotation_pie.py` — generic HOMER annotation composition pies
