#!/usr/bin/env python3
"""Scale each bigWig sample column-block in a computeMatrix .gz independently (deeptools format).

Copied from CUT&RUN `analysis/p65_ER_cobinding/signal_profiles/normalize_matrix_samples.py`.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from deeptools.heatmapper import heatmapper


def _denominator(block: np.ma.MaskedArray, mode: str) -> float | None:
    vals = block.compressed()
    if vals.size == 0:
        return None
    if mode == "mean":
        d = float(np.mean(vals))
    elif mode == "max":
        d = float(np.max(vals))
    elif mode == "p99":
        d = float(np.percentile(vals, 99))
    elif mode == "p95":
        d = float(np.percentile(vals, 95))
    else:
        raise ValueError(f"Unknown mode: {mode}")
    if not np.isfinite(d) or d <= 0:
        return None
    return d


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--input", type=str, required=True)
    p.add_argument("-o", "--output", type=str, required=True)
    p.add_argument(
        "--mode",
        choices=("none", "mean", "max", "p95", "p99"),
        default="p99",
        help="Divide each sample's columns by this statistic (computed over all regions/bins).",
    )
    args = p.parse_args()

    hm = heatmapper()
    hm.read_matrix_file(args.input)
    m = hm.matrix.matrix
    bounds = hm.matrix.sample_boundaries

    for si in range(len(bounds) - 1):
        lo, hi = bounds[si], bounds[si + 1]
        block = m[:, lo:hi]
        if args.mode == "none":
            continue
        den = _denominator(block, args.mode)
        if den is None:
            continue
        m[:, lo:hi] /= den

    hm.save_matrix(args.output)
    print(f"Wrote {args.output} (mode={args.mode})", file=sys.stderr)


if __name__ == "__main__":
    main()

