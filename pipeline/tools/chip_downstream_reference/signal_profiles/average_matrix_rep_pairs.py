#!/usr/bin/env python3
"""Average replicate sample blocks in a deepTools computeMatrix (.gz) into condition-level blocks.

This is used as a robust alternative to bigwigCompare --operation mean when bigwigCompare is slow/hangs.

Example pairs:
  Veh,Veh_rep1,Veh_rep2
  E2,E2_rep1,E2_rep2
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable

import numpy as np
from deeptools.heatmapper import heatmapper


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--input", required=True, help="computeMatrix output .gz")
    p.add_argument("-o", "--output", required=True, help="averaged matrix .gz")
    p.add_argument(
        "--pairs",
        nargs="+",
        required=True,
        help="Triples: OUT_LABEL,REP1_LABEL,REP2_LABEL (comma-separated)",
    )
    return p.parse_args()


def _triples(pairs: Iterable[str]) -> list[tuple[str, str, str]]:
    out: list[tuple[str, str, str]] = []
    for s in pairs:
        parts = [x.strip() for x in s.split(",")]
        if len(parts) != 3 or not all(parts):
            raise ValueError(f"Bad --pairs entry (need OUT,REP1,REP2): {s!r}")
        out.append((parts[0], parts[1], parts[2]))
    return out


def main() -> None:
    args = parse_args()
    triples = _triples(args.pairs)

    hm = heatmapper()
    hm.read_matrix_file(args.input)

    m = hm.matrix.matrix  # masked array
    bounds = list(hm.matrix.sample_boundaries)
    labels = list(getattr(hm.matrix, "sample_labels", []))
    if not labels:
        print("Matrix missing sample_labels; cannot average by label.", file=sys.stderr)
        raise SystemExit(1)

    label_to_idx = {lab: i for i, lab in enumerate(labels)}

    blocks: list[np.ma.MaskedArray] = []
    new_labels: list[str] = []
    for out_lab, rep1, rep2 in triples:
        if rep1 not in label_to_idx or rep2 not in label_to_idx:
            print(f"Missing rep labels in matrix: {rep1!r} or {rep2!r}", file=sys.stderr)
            raise SystemExit(1)
        i1 = label_to_idx[rep1]
        i2 = label_to_idx[rep2]
        lo1, hi1 = bounds[i1], bounds[i1 + 1]
        lo2, hi2 = bounds[i2], bounds[i2 + 1]
        if (hi1 - lo1) != (hi2 - lo2):
            print(f"Rep blocks have different widths: {rep1} vs {rep2}", file=sys.stderr)
            raise SystemExit(1)
        b1 = m[:, lo1:hi1]
        b2 = m[:, lo2:hi2]
        blocks.append((b1 + b2) / 2.0)
        new_labels.append(out_lab)

    new_m = np.ma.concatenate(blocks, axis=1)

    # Update heatmapper matrix in-place
    hm.matrix.matrix = new_m
    step = blocks[0].shape[1]
    hm.matrix.sample_boundaries = [i * step for i in range(len(blocks) + 1)]
    hm.matrix.sample_labels = new_labels

    hm.save_matrix(args.output)
    print(f"Wrote {args.output} (averaged {len(triples)} pairs)", file=sys.stderr)


if __name__ == "__main__":
    main()

