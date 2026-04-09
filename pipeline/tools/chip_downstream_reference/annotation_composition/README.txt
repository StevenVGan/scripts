Method 2 peak genomic annotation composition
===========================================

Purpose
-------
Summarize HOMER-style genomic annotations for:
  - Baseline (Arm A, single-ligand consensus): p65 TNF + ER E2
  - Method 2 cobinding peaks: arms A (cross), B (cotreat), C (cotreat-exclusive)

Cobinding tables are p65-centric (annotatePeaks on p65 intervals). Bar chart legend uses
Cross-condition, Cotreatment, Cotreat-exclusive (same spirit as conservation plots), plus
p65 consensus (TNF) and ER consensus (E2). ER in the bar plot is a crude reference (different
feature than p65 cobinding). All PNGs for a caller: output/figures/<caller>/.

Inputs
------
Under analysis/p65_ER_cobinding/<macs|homer>/:
  method2_A_cross/p65_TNF_consensus.annotatePeaks.txt
  method2_A_cross/ER_E2_consensus.bed  +  peaks/<macs3|homer>/ER_E2_rep1.annotatePeaks.txt
  method2_{A_cross,B_cotreat,C_cotreat_exclusive}/method2_all_p65_nearestER_within2kb.tsv

Arm C TSV must include the annotatePeaks header (Phase 0 fix in run_method2_cobinding.sh).

Run
---
  cd analysis/p65_ER_cobinding/annotation_composition
  ./run_annotation_composition.sh

Or (from this directory):
  python3 method2_annotation_composition.py --cobinding-root .. --chip-root ../../..
  # chip-root = GSE59530 project root (parent of analysis/), three levels up from here.

Options: --caller macs|homer|both  --no-er-baseline  --out-dir ./output  --dpi 200

Outputs
-------
See OUTPUT_LAYOUT.txt and output/README.txt. All products are under output/ by default:
  tables/<caller>/counts_long.tsv, counts_wide_p65.tsv
  figures/<caller>/pie_*.png, bars_baseline_vs_cobinding_ABC.png
  run_manifest.tsv
