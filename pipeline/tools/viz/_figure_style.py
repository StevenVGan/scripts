"""Publication-quality matplotlib styling — shared lab default.

Sized for 14-22-inch canvases so labels stay legible when figures are reduced to
print scale; vector-editable text in PDFs (pdf.fonttype=42), sans-serif. Promoted
from the Ren/Reed multiome analyses (the "big-canvas" variant).

Usage:
    import matplotlib
    matplotlib.use("Agg")
    from _figure_style import apply_publication_style
    apply_publication_style()
    import matplotlib.pyplot as plt
"""
from __future__ import annotations

import matplotlib as mpl


def apply_publication_style() -> None:
    mpl.rcParams.update({
        "font.size":        16,
        "axes.titlesize":   20,
        "axes.labelsize":   17,
        "xtick.labelsize":  14,
        "ytick.labelsize":  14,
        "legend.fontsize":  14,
        "figure.titlesize": 22,
        "font.family":      "sans-serif",
        "font.sans-serif":  ["DejaVu Sans", "Arial", "Helvetica"],
        "pdf.fonttype":     42,
        "ps.fonttype":      42,
        "svg.fonttype":     "none",
        "savefig.dpi":      300,
        "figure.dpi":       110,
        "savefig.bbox":     "tight",
        "axes.linewidth":   0.9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.direction":  "out",
        "ytick.direction":  "out",
        "legend.frameon":   False,
    })
