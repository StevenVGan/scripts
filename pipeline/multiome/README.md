# Multiome (10x snATAC + snRNA) pipeline

**Status:** intentionally empty.

The first multiome project (`seq/multiome/Multiome_GSE278576_AgingHippocampus_Ren/`)
is being built with project-local Python scripts under its own `script/local/`.
No shared pipeline lives here yet because we have only one project — promoting
prematurely would design a generic interface against a single dataset.

Per CONVENTIONS §5, scripts are promoted into `scripts/pipeline/<assay>/` (or
`scripts/pipeline/tools/`) once the same script is needed in a 2nd project.
When a 2nd multiome project arrives, expect candidates like:

- a generic 10x cellranger-arc H5 → AnnData loader
- per-cell QC (n_genes, %mito, n_fragments, TSS enrichment, FRiP)
- doublet calling (Scrublet + AMULET)
- donor-batch integration (harmony / scVI)
- SnapATAC2 spectral + clustering wrapper
- joint analysis (muon WNN, peak-gene linkage)

Until then: see `seq/multiome/Multiome_GSE278576_AgingHippocampus_Ren/script/local/`
for the working code.

## Env

`scripts/env/sc.yml` (separate from `bio` to avoid bulk-seq vs. SC dep
conflicts; first SC project in the workspace).
