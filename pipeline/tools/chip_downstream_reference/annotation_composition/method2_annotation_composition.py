#!/usr/bin/env python3
"""Genomic annotation composition: Method 2 cobinding A/B/C vs Arm A baselines (p65 TNF, ER E2).

Reads HOMER/MACS annotatePeaks-style tables, buckets Annotation strings, writes long/wide TSVs
and figures under output/ (see OUTPUT_LAYOUT.txt).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Ordered categories for consistent colors and bar order
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

# Distinct colors (tab20-based, enough for 9 categories)
_COLORS = plt.cm.tab20(np.linspace(0, 1, 20))
CATEGORY_COLORS = {cat: _COLORS[i % 20] for i, cat in enumerate(CATEGORY_ORDER)}

# Pie titles / bar legend (match conservation-style arm wording; see downstream/conservation REGION_PLOT_LABELS)
DISPLAY_LABEL = {
    "p65_consensus_A": "p65 consensus (TNF)",
    "ER_consensus_A": "ER consensus (E2)",
    "cobinding_A": "Cross-condition",
    "cobinding_B": "Cotreatment",
    "cobinding_C": "Cotreat-exclusive",
}

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
    with path.open(newline="") as f:
        return list(csv.reader(f, delimiter="\t"))


def annotations_from_annotate_peaks(path: Path) -> list[str]:
    rows = _read_tsv_rows(path)
    if not rows:
        return []
    header = rows[0]
    if "Annotation" not in header:
        raise ValueError(f"No Annotation column in {path}")
    idx = header.index("Annotation")
    out: list[str] = []
    for row in rows[1:]:
        if len(row) > idx:
            out.append(row[idx])
    return out


def annotations_from_cobinding_tsv(path: Path) -> list[str]:
    rows = _read_tsv_rows(path)
    if not rows:
        return []
    if "Annotation" in rows[0]:
        idx = rows[0].index("Annotation")
        return [r[idx] for r in rows[1:] if len(r) > idx]
    idx = 7
    anns: list[str] = []
    for row in rows:
        if len(row) > idx:
            anns.append(row[idx])
    return anns


def load_er_consensus_peak_ids(bed_path: Path) -> set[str]:
    ids: set[str] = set()
    for line in bed_path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            ids.add(parts[3].strip())
    return ids


def annotations_er_consensus(bed_path: Path, annotate_path: Path) -> list[str]:
    peak_ids = load_er_consensus_peak_ids(bed_path)
    rows = _read_tsv_rows(annotate_path)
    if not rows:
        raise ValueError(f"Empty annotate: {annotate_path}")
    h = rows[0]
    i_peak = 0  # HOMER/MACS header col0 is PeakID (often "PeakID (cmd=...)")
    try:
        i_ann = h.index("Annotation")
    except ValueError:
        i_ann = 8
    out: list[str] = []
    for row in rows[1:]:
        if len(row) <= max(i_peak, i_ann):
            continue
        pid = row[i_peak].strip()
        if pid in peak_ids:
            out.append(row[i_ann])
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
    """Categories in CATEGORY_ORDER first, then any extras."""
    seen = set()
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


def load_expected_peaks(cobinding_root: Path) -> dict[tuple[str, str], int]:
    """Map (caller_lower, arm_slug) -> n_peaks from summary/method2_ABC_counts.tsv."""
    summary = cobinding_root / "summary" / "method2_ABC_counts.tsv"
    expected: dict[tuple[str, str], int] = {}
    if not summary.is_file():
        return expected
    with summary.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            caller = row["caller"].strip()
            arm = row["arm"].strip()
            np_ = int(row["n_peaks"])
            cl = "macs" if caller == "MACS3" else caller.lower()
            expected[(cl, arm)] = np_
    return expected


COBINDING_ARM_KEY = {
    "cobinding_A": "A_cross",
    "cobinding_B": "B_cotreat",
    "cobinding_C": "C_cotreat_exclusive",
}


def validate_cobinding_counts(
    caller: str, set_id: str, n: int, expected: dict[tuple[str, str], int]
) -> None:
    if not set_id.startswith("cobinding_"):
        return
    arm = COBINDING_ARM_KEY.get(set_id)
    if not arm:
        return
    exp = expected.get((caller, arm))
    if exp is None:
        return
    if n != exp:
        print(
            f"[WARN] {caller} {set_id}: parsed n={n} but summary/method2_ABC_counts.tsv "
            f"expects n_peaks={exp}",
            file=sys.stderr,
        )


def plot_pie(labels: list[str], values: list[int], title: str, out_path: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    colors = [CATEGORY_COLORS.get(l, (0.5, 0.5, 0.5, 1.0)) for l in labels]
    total = sum(values)
    if total == 0:
        ax.text(0.5, 0.5, "n=0", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        wedges, texts, autotexts = ax.pie(
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
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_bars_annotation_comparison(
    series: dict[str, Counter],
    out_path: Path,
    dpi: int,
    merge_under_pct: float,
    include_er_in_bars: bool,
) -> None:
    """Grouped bars: one group per category; series includes baselines + cobinding (optional ER)."""
    order: list[str] = []
    if "p65_consensus_A" in series:
        order.append("p65_consensus_A")
    if include_er_in_bars and "ER_consensus_A" in series:
        order.append("ER_consensus_A")
    for k in ["cobinding_A", "cobinding_B", "cobinding_C"]:
        if k in series:
            order.append(k)
    if len(order) < 2:
        return
    merged_series = {k: merge_small_categories(series[k], merge_under_pct) for k in order}
    all_cats: set[str] = set()
    for c in merged_series.values():
        all_cats.update(c.keys())
    cat_list = [c for c in CATEGORY_ORDER if c in all_cats]
    cat_list += sorted(all_cats - set(cat_list))

    n_bars = len(order)
    x = np.arange(len(cat_list))
    width = min(0.22, 0.9 / max(n_bars, 1))
    fig_w = max(11, len(cat_list) * (0.45 + 0.12 * n_bars))
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    for i, key in enumerate(order):
        counts = merged_series[key]
        fracs = [counts.get(c, 0) / sum(counts.values()) if sum(counts.values()) else 0 for c in cat_list]
        offset = (i - (n_bars - 1) / 2) * width
        lab = DISPLAY_LABEL.get(key, key)
        ax.bar(x + offset, fracs, width, label=lab, color=plt.cm.Set2(i / max(n_bars, 1)))

    ax.set_ylabel("Fraction of peaks")
    ax.set_title(
        "Genomic annotation: p65 / ER baselines vs Method 2 cobinding\n"
        "(ER is a crude reference; cobinding annotations are p65-centric)"
    )
    ax.set_xticks(x)
    ax.set_xticklabels(cat_list, rotation=35, ha="right", fontsize=8)
    ax.legend(loc="upper right", fontsize=7)
    ax.set_ylim(0, min(1.0, ax.get_ylim()[1]))
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_counts_long(
    rows: list[dict],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["caller", "set_id", "feature", "category", "count", "fraction", "n_peaks"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in fields})


def write_counts_wide_p65(rows_long: list[dict], path: Path) -> None:
    p65_sets = {"p65_consensus_A", "cobinding_A", "cobinding_B", "cobinding_C"}
    by_set: dict[str, dict[str, float]] = {}
    for row in rows_long:
        if row["feature"] != "p65" or row["set_id"] not in p65_sets:
            continue
        sid = row["set_id"]
        by_set.setdefault(sid, {})[row["category"]] = float(row["fraction"])

    all_cats = set()
    for d in by_set.values():
        all_cats.update(d)
    cat_list = [c for c in CATEGORY_ORDER if c in all_cats] + sorted(all_cats - set(CATEGORY_ORDER))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["set_id"] + cat_list)
        for sid in ["p65_consensus_A", "cobinding_A", "cobinding_B", "cobinding_C"]:
            if sid not in by_set:
                continue
            w.writerow([sid] + [f"{by_set[sid].get(c, 0):.6f}" for c in cat_list])


def append_manifest(out_root: Path, argv: list[str], cobinding_root: Path, chip_root: Path) -> None:
    man = out_root / "run_manifest.tsv"
    row = {
        "iso_timestamp": datetime.now(timezone.utc).isoformat(),
        "cobinding_root": str(cobinding_root),
        "chip_root": str(chip_root),
        "argv": " ".join(argv),
    }
    write_header = not man.is_file()
    with man.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()), delimiter="\t", lineterminator="\n")
        if write_header:
            w.writeheader()
        w.writerow(row)


def run_caller(
    caller: str,
    cobinding_root: Path,
    chip_root: Path,
    out_root: Path,
    peaks_subdir: str,
    include_er: bool,
    merge_under_pct: float,
    dpi: int,
    expected: dict[tuple[str, str], int],
) -> list[dict]:
    cob = cobinding_root / caller
    peaks_dir = chip_root / "peaks" / peaks_subdir

    configs: list[tuple[str, str, Path | None, Path | None, str]] = [
        (
            "p65_consensus_A",
            "p65",
            cob / "method2_A_cross" / "p65_TNF_consensus.annotatePeaks.txt",
            None,
            "annotate",
        ),
        (
            "cobinding_A",
            "p65",
            cob / "method2_A_cross" / "method2_all_p65_nearestER_within2kb.tsv",
            None,
            "cobinding_tsv",
        ),
        (
            "cobinding_B",
            "p65",
            cob / "method2_B_cotreat" / "method2_all_p65_nearestER_within2kb.tsv",
            None,
            "cobinding_tsv",
        ),
        (
            "cobinding_C",
            "p65",
            cob / "method2_C_cotreat_exclusive" / "method2_all_p65_nearestER_within2kb.tsv",
            None,
            "cobinding_tsv",
        ),
    ]
    if include_er:
        configs.append(
            (
                "ER_consensus_A",
                "ER",
                cob / "method2_A_cross" / "ER_E2_consensus.bed",
                peaks_dir / "ER_E2_rep1.annotatePeaks.txt",
                "er_bed",
            )
        )

    long_rows: list[dict] = []
    series_for_bars: dict[str, Counter] = {}

    fig_dir = out_root / "figures" / caller
    fig_dir.mkdir(parents=True, exist_ok=True)

    for set_id, feature, p1, p2, kind in configs:
        if kind == "annotate":
            assert p1 is not None
            anns = annotations_from_annotate_peaks(p1)
        elif kind == "cobinding_tsv":
            assert p1 is not None
            anns = annotations_from_cobinding_tsv(p1)
            validate_cobinding_counts(caller, set_id, len(anns), expected)
        elif kind == "er_bed":
            assert p1 is not None and p2 is not None
            anns = annotations_er_consensus(p1, p2)
        else:
            raise ValueError(kind)

        raw_counts = Counter(classify(a) for a in anns)
        counts = merge_small_categories(raw_counts, merge_under_pct)
        n_peaks = len(anns)

        for cat, c in counts.items():
            frac = c / n_peaks if n_peaks else 0.0
            long_rows.append(
                {
                    "caller": caller,
                    "set_id": set_id,
                    "feature": feature,
                    "category": cat,
                    "count": c,
                    "fraction": f"{frac:.6f}",
                    "n_peaks": n_peaks,
                }
            )

        labels, values = counter_to_ordered(counts)
        disp = DISPLAY_LABEL.get(set_id, set_id)
        title = f"{caller}\n{disp}\nn={n_peaks}"
        plot_pie(labels, values, title, fig_dir / f"pie_{set_id}.png", dpi)
        if set_id in ("p65_consensus_A", "ER_consensus_A", "cobinding_A", "cobinding_B", "cobinding_C"):
            series_for_bars[set_id] = counts

    plot_bars_annotation_comparison(
        series_for_bars,
        fig_dir / "bars_baseline_vs_cobinding_ABC.png",
        dpi,
        merge_under_pct,
        include_er_in_bars=include_er,
    )

    tab_dir = out_root / "tables" / caller
    write_counts_long(long_rows, tab_dir / "counts_long.tsv")
    write_counts_wide_p65(long_rows, tab_dir / "counts_wide_p65.tsv")

    return long_rows


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    default_root = here.parent  # p65_ER_cobinding
    default_chip = default_root.parent.parent  # .../MCF7_ER_p65_ChIP_GSE59530 (has peaks/)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cobinding-root", type=Path, default=default_root)
    p.add_argument("--chip-root", type=Path, default=default_chip)
    p.add_argument("--out-dir", type=Path, default=here / "output")
    p.add_argument("--caller", choices=("macs", "homer", "both"), default="both")
    p.add_argument(
        "--no-er-baseline",
        action="store_true",
        help="Skip ER consensus pies/tables and omit ER from the bar chart.",
    )
    p.add_argument("--merge-small-pct", type=float, default=0.0, help="Merge categories below this %% into Other (pies/bars).")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cobinding_root = args.cobinding_root.resolve()
    chip_root = args.chip_root.resolve()
    out_root = args.out_dir.resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    append_manifest(out_root, sys.argv, cobinding_root, chip_root)
    expected = load_expected_peaks(cobinding_root)

    callers = ["macs", "homer"] if args.caller == "both" else [args.caller]
    peaks_sub = {"macs": "macs3", "homer": "homer"}

    for caller in callers:
        run_caller(
            caller,
            cobinding_root,
            chip_root,
            out_root,
            peaks_sub[caller],
            include_er=not args.no_er_baseline,
            merge_under_pct=args.merge_small_pct,
            dpi=args.dpi,
            expected=expected,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
