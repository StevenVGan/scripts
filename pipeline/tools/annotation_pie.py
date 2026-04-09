#!/usr/bin/env python3
"""Pie chart of genomic annotation composition (HOMER-style Annotation strings).

Input formats:
  homer  — tab-delimited table with an 'Annotation' header column (annotatePeaks.txt)
  tsv    — tab-delimited; specify column by name (--annotation-column) or 1-based index
  lines  — one raw annotation string per line (same text HOMER would put in Annotation)

Categories follow the same bucketing as the GSE59530 Method 2 composition script
(promoter / intron / intergenic / …).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CATEGORY_ORDER = [
    "Promoter (incl. TSS)",
    "Intronic",
    "Intergenic",
    "Exonic",
    "TTS (transcription end)",
    "Non-coding annotation",
    "Repeat / satellite / other feature",
    "NA / unannotated",
    "Other",
]

_COLORS = plt.cm.tab20(np.linspace(0, 1, 20))
CATEGORY_COLORS = {cat: _COLORS[i % 20] for i, cat in enumerate(CATEGORY_ORDER)}


def classify(ann: str) -> str:
    if not ann or ann.strip() in ("", "NA"):
        return "NA / unannotated"
    a = ann.strip()
    low = a.lower()
    if low.startswith("promoter"):
        return "Promoter (incl. TSS)"
    if low.startswith("intron"):
        return "Intronic"
    if low.startswith("exon"):
        return "Exonic"
    if low.startswith("tts"):
        return "TTS (transcription end)"
    if low.startswith("intergenic"):
        return "Intergenic"
    if "non-coding" in low or low.startswith("noncoding"):
        return "Non-coding annotation"
    if "|" in a or "satellite" in low or "simple" in low:
        return "Repeat / satellite / other feature"
    return "Other"


def _read_tsv_rows(path: Path) -> list[list[str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.reader(f, delimiter="\t"))


def _header_index_ci(header: list[str], name: str) -> int:
    target = name.strip().lower()
    for i, h in enumerate(header):
        if h.strip().lower() == target:
            return i
    raise ValueError(f"No {name!r} column in header of TSV (case-insensitive match)")


def load_annotations_homer(path: Path) -> list[str]:
    rows = _read_tsv_rows(path)
    if not rows:
        return []
    header = rows[0]
    idx = _header_index_ci(header, "Annotation")
    out: list[str] = []
    for row in rows[1:]:
        if len(row) > idx:
            out.append(row[idx])
    return out


def load_annotations_tsv_column(path: Path, column: str) -> list[str]:
    rows = _read_tsv_rows(path)
    if not rows:
        return []
    header = rows[0]
    if column.isdigit():
        idx = int(column) - 1
        if idx < 0 or idx >= len(header):
            raise ValueError(f"Column index {column} out of range (n_cols={len(header)})")
    else:
        lower_map = {h.lower(): i for i, h in enumerate(header)}
        key = column.lower()
        if key not in lower_map:
            raise ValueError(f"Column {column!r} not in header: {header[:15]}...")
        idx = lower_map[key]
    out: list[str] = []
    for row in rows[1:]:
        if len(row) > idx:
            out.append(row[idx])
    return out


def load_annotations_lines(path: Path) -> list[str]:
    out: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            out.append(s)
    return out


def merge_small_categories(counts: Counter, merge_under_pct: float) -> Counter:
    if merge_under_pct <= 0:
        return counts
    total = sum(counts.values())
    if total == 0:
        return counts
    merged = Counter()
    for cat, c in counts.items():
        pct = 100.0 * c / total
        if pct < merge_under_pct and cat != "Other":
            merged["Other"] += c
        else:
            merged[cat] += c
    return merged


def counter_to_ordered(counts: Counter) -> tuple[list[str], list[int]]:
    seen: set[str] = set()
    labels: list[str] = []
    values: list[int] = []
    for cat in CATEGORY_ORDER:
        if cat in counts:
            labels.append(cat)
            values.append(counts[cat])
            seen.add(cat)
    for cat, c in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        if cat not in seen:
            labels.append(cat)
            values.append(c)
    return labels, values


def plot_pie(labels: list[str], values: list[int], title: str, out_path: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    colors = [CATEGORY_COLORS.get(l, (0.5, 0.5, 0.5, 1.0)) for l in labels]
    total = sum(values)
    if total == 0:
        ax.text(0.5, 0.5, "n=0", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        wedges, _texts, _autotexts = ax.pie(
            values,
            labels=None,
            autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
            colors=colors,
            pctdistance=0.75,
        )
        ax.legend(
            wedges,
            [f"{l} ({100 * v / total:.1f}%)" for l, v in zip(labels, values)],
            title="Category",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            fontsize=8,
        )
    ax.set_title(title)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def detect_format(path: Path) -> str:
    """Heuristic: if first line is TSV and contains Annotation header -> homer."""
    try:
        with path.open(encoding="utf-8", errors="replace") as f:
            line = f.readline()
    except OSError:
        return "lines"
    if "\t" not in line:
        return "lines"
    fields = line.rstrip("\n").split("\t")
    if any(f.strip().lower() == "annotation" for f in fields):
        return "homer"
    return "tsv"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", type=Path, required=True)
    p.add_argument(
        "--format",
        choices=("auto", "homer", "tsv", "lines"),
        default="auto",
        help="auto: Annotation column in header -> homer; else tab file -> need --annotation-column; no tabs -> lines",
    )
    p.add_argument(
        "--annotation-column",
        default="Annotation",
        help='For --format tsv (and auto when not homer): column name or 1-based index (default "Annotation")',
    )
    p.add_argument("--out", "-o", type=Path, required=True, help="Output PNG path")
    p.add_argument("--title", default="", help="Plot title (default: input basename)")
    p.add_argument("--merge-small-pct", type=float, default=0.0, help="Merge categories below this %% into Other")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument(
        "--counts-tsv",
        type=Path,
        default=None,
        help="Optional path to write category counts (tab-separated)",
    )
    args = p.parse_args()

    inp = args.input.resolve()
    if not inp.is_file():
        print(f"ERROR: not a file: {inp}", file=sys.stderr)
        return 1

    fmt = args.format
    if fmt == "auto":
        fmt = detect_format(inp)
        if fmt == "tsv" and args.annotation_column.lower() == "annotation":
            # Tab file but no Annotation in header — user must pass real column
            rows = _read_tsv_rows(inp)
            if rows:
                try:
                    _header_index_ci(rows[0], "Annotation")
                except ValueError:
                    print(
                        "[WARN] auto detected 'tsv' but no Annotation column; "
                        "using --annotation-column as given (may fail if invalid).",
                        file=sys.stderr,
                    )

    try:
        if fmt == "homer":
            raw = load_annotations_homer(inp)
        elif fmt == "tsv":
            raw = load_annotations_tsv_column(inp, args.annotation_column)
        else:
            raw = load_annotations_lines(inp)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    counts: Counter = Counter()
    for a in raw:
        counts[classify(a)] += 1

    counts = merge_small_categories(counts, args.merge_small_pct)
    labels, values = counter_to_ordered(counts)

    title = args.title.strip() or f"{inp.name} (n={sum(values)})"
    plot_pie(labels, values, title, args.out.resolve(), args.dpi)
    print(f"Wrote {args.out} (n_peaks={sum(values)})", file=sys.stderr)

    if args.counts_tsv:
        outp = args.counts_tsv.resolve()
        outp.parent.mkdir(parents=True, exist_ok=True)
        with outp.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["category", "count"])
            for lab, val in zip(labels, values):
                w.writerow([lab, val])
        print(f"Wrote {outp}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
