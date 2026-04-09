ER vs p65 ChIP-seq signal profiles / heatmaps (HOMER peaks, 4 conditions, hg38)
==============================================================================

Purpose
-------
Mirror the CUT&RUN `signal_profiles/` workflow, but for GSE59530 ChIP-seq and 4 conditions
(Veh, E2, TNF, E2TNF). We plot average bigWig signal around three region classes using
deepTools `computeMatrix` + `plotHeatmap`. Outputs are separate figures
`figure_ER_by_siteclass.png` and `figure_p65_by_siteclass.png`. ER signal is much
stronger than p65, so `run_signal_profiles_chip.sh` sets a higher heatmap/profile
ceiling for ER than for p65 (`--zMax` / `--yMax`, defaults 75 vs 15). Override with
environment variables `ZMAX_ER`, `YMAX_ER`, `ZMAX_P65`, `YMAX_P65` if needed.

Region classes (HOMER-based)
----------------------------
- Cotreat_cobind: Method2 cotreatment cobinding (Arm B; HOMER), derived from
  `../downstream/lists/homer/method2_B_cotreat.bed` (heatmap row label)
- p65_only: p65_TNF consensus peaks (HOMER) minus Cotreat_cobind
- ER_only_10pct: ER_E2 consensus peaks (HOMER) with no overlap of p65_TNF consensus expanded by PROX,
  then reproducible 10% subsample for plotting

Tracks (4 condition columns)
----------------------------
Input bigWigs (replicate tracks) live under:
  `../../../align/track/`

Replicate bigWigs are averaged inside the matrix (`average_matrix_rep_pairs.py`) after `computeMatrix`.

Run
---
  conda activate bio    # bedtools, deepTools, python3
  bash run_signal_profiles_chip.sh

  # If deepTools is only installed in the env and not on PATH:
  # conda run -n bio bash run_signal_profiles_chip.sh

Outputs (default under this directory)
--------------------------------------
  tracks_manifest.tsv
  tracks_mean/*.bw
  cobind_E2TNF.bed, p65_TNF_only.bed, ER_E2_only_10pct.bed
  site_bed_qc.txt

  matrix_ER_by_siteclass.pre_norm.gz
  matrix_ER_by_siteclass.gz
  figure_ER_by_siteclass.png

  matrix_p65_by_siteclass.pre_norm.gz
  matrix_p65_by_siteclass.gz
  figure_p65_by_siteclass.png

Normalization
-------------
`normalize_matrix_samples.py` scales each sample independently (default: p99) so different
conditions are on comparable color scales.

