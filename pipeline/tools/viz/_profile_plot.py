"""_profile_plot.py — custom matplotlib profile plotter for deepTools matrices.

Replaces deepTools `plotProfile` so we can:
  1. Place the legend in the **rightmost** subplot (deepTools always puts it
     in the first/leftmost subplot regardless of --legendLocation).
  2. Honor the shared `_figure_style.py` rcParams (deepTools subprocesses
     don't pick those up).
  3. Render exactly the --perGroup layout (one subplot per region group,
     samples overlaid).

Matrix format (deepTools computeMatrix output, gzipped tab-separated):
  - Line 1: `@<json>` header containing:
      group_labels, group_boundaries (row indices delimiting region groups)
      sample_labels, sample_boundaries (column indices delimiting samples)
      upstream, downstream, bin size (per-sample arrays)
  - Lines 2+: chrom, start, end, name, score, strand, then matrix values

Usage:
    from _profile_plot import plot_profile
    plot_profile(mat_gz, region_labels, sample_labels,
                 out_prefix=Path("figures/.../atac_signal_profile"),
                 plot_height=5, plot_width=5)
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _read_matrix(mat_gz: str | Path):
    with gzip.open(str(mat_gz), "rt") as f:
        header_line = f.readline()
        assert header_line.startswith("@"), f"unexpected header: {header_line[:40]}"
        header = json.loads(header_line[1:])
        rows = []
        for line in f:
            cols = line.rstrip("\n").split("\t")
            # cols 0..5 = chrom, start, end, name, score, strand; values from col 6.
            vals = [np.nan if v in ("nan", "NA", "") else float(v) for v in cols[6:]]
            rows.append(vals)
    return header, np.array(rows, dtype=np.float64)


def plot_profile(mat_gz, region_labels, sample_labels, out_prefix,
                 plot_height: float = 5.0, plot_width: float = 5.0,
                 sample_colors=None, ref_label: str = "center"):
    """Render --perGroup profile with legend on the rightmost subplot.

    region_labels : list[str], one per region group (BED) in computeMatrix.
    sample_labels : list[str], one per sample (bigwig) in computeMatrix.
    out_prefix    : Path or str, output prefix (extension appended).
    """
    header, data = _read_matrix(mat_gz)
    sample_boundaries = header["sample_boundaries"]
    group_boundaries = header["group_boundaries"]
    n_samples = len(sample_labels)
    n_groups = len(region_labels)
    bin_size = header["bin size"][0]
    upstream = header["upstream"][0]
    downstream = header["downstream"][0]
    n_bins = sample_boundaries[1] - sample_boundaries[0]
    # x axis: position relative to ref_point in bp, centered on 0
    x = np.arange(-upstream, downstream, bin_size) + bin_size / 2

    fig, axes = plt.subplots(
        1, n_groups, figsize=(plot_width * n_groups, plot_height + 0.5),
        sharey=True,
    )
    if n_groups == 1:
        axes = [axes]

    if sample_colors is None:
        cmap = plt.get_cmap("tab10")
        sample_colors = [cmap(i % 10) for i in range(n_samples)]

    for gi, ax in enumerate(axes):
        r0, r1 = group_boundaries[gi], group_boundaries[gi + 1]
        for si in range(n_samples):
            c0, c1 = sample_boundaries[si], sample_boundaries[si + 1]
            slab = data[r0:r1, c0:c1]
            mean = np.nanmean(slab, axis=0)
            # Guard against shape mismatch if the file has slightly different bin count
            mx = min(len(x), len(mean))
            ax.plot(x[:mx], mean[:mx], color=sample_colors[si], lw=2.2,
                    label=sample_labels[si])
        ax.set_title(f"{region_labels[gi]}  (n={r1 - r0:,})")
        ax.set_xlabel(f"Distance from {ref_label} (bp)")
        if gi == 0:
            ax.set_ylabel("Mean signal")
        ax.axvline(0, color="grey", lw=0.6, ls="--", alpha=0.6)
        # Legend ONLY in the rightmost subplot
        if gi == n_groups - 1:
            ax.legend(loc="upper right", frameon=True, framealpha=0.9,
                      edgecolor="none")

    fig.tight_layout()
    out_prefix = str(out_prefix)
    out_paths = []
    for ext in ("png", "pdf"):
        out = f"{out_prefix}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        out_paths.append(out)
    plt.close(fig)
    return out_paths
