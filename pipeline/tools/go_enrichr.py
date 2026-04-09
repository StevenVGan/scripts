#!/usr/bin/env python3
"""Enrichr enrichment + horizontal bar plot of top terms.

Input modes:
  genes   — one gene symbol per line (optional # comments, blank lines skipped)
  annotate — HOMER/MACS annotatePeaks-style TSV: header with Gene Name (or --gene-column)
  bed     — BED with ≥4 columns; column 4 (name) is treated as gene symbol per row

Requires: gseapy, matplotlib, and network access to the Enrichr web API.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import tempfile
import time
from pathlib import Path

DEFAULT_GENE_SET = "GO_Biological_Process_2023"


def load_genes_plain(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        names.append(s.upper())
    return sorted(set(names))


def _read_tsv(path: Path) -> list[list[str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.reader(f, delimiter="\t"))


def _split_gene_cell(cell: str) -> list[str]:
    """Split 'GENE1, GENE2' / 'GENE1; GENE2' into symbols."""
    parts = re.split(r"[,;/]\s*", cell.strip())
    return [p.strip().upper() for p in parts if p.strip()]


def load_genes_annotate(path: Path, gene_column: str) -> list[str]:
    rows = _read_tsv(path)
    if not rows:
        return []
    header = rows[0]
    col = gene_column.strip()
    try:
        if col.isdigit():
            idx = int(col) - 1
            if idx < 0:
                raise ValueError
        else:
            # case-insensitive match
            lower_map = {h.lower(): i for i, h in enumerate(header)}
            key = col.lower()
            if key not in lower_map:
                raise KeyError(col)
            idx = lower_map[key]
    except (KeyError, ValueError) as e:
        raise SystemExit(
            f"Column {gene_column!r} not found in header of {path}. "
            f"Columns: {header[:20]}{'...' if len(header) > 20 else ''}"
        ) from e

    names: list[str] = []
    for row in rows[1:]:
        if len(row) <= idx:
            continue
        cell = row[idx].strip()
        if not cell or cell == "NA":
            continue
        names.extend(_split_gene_cell(cell))
    return sorted(set(names))


def load_genes_bed(path: Path, gene_col: int) -> list[str]:
    """gene_col is 1-based BED column index (default 4 = name)."""
    idx = gene_col - 1
    names: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("track ") or s.startswith("browser "):
            continue
        parts = s.split("\t")
        if len(parts) <= idx:
            continue
        g = parts[idx].strip()
        if not g or g == ".":
            continue
        names.extend(_split_gene_cell(g))
    return sorted(set(names))


def run_enrichr(genes: list[str], gene_set: str, organism: str, tmp_parent: Path):
    import gseapy as gp

    tmpd = tmp_parent / "enr"
    tmpd.mkdir(parents=True, exist_ok=True)
    return gp.enrichr(
        gene_list=genes,
        gene_sets=[gene_set],
        organism=organism,
        no_plot=True,
        outdir=str(tmpd),
        verbose=False,
    )


def run_enrichr_with_retry(
    genes: list[str],
    gene_set: str,
    organism: str,
    tmp_parent: Path,
    retries: int,
    retry_wait_sec: int,
):
    last_exc = None
    for i in range(retries + 1):
        try:
            return run_enrichr(genes, gene_set, organism, tmp_parent)
        except Exception as e:
            last_exc = e
            if i >= retries:
                raise
            time.sleep(retry_wait_sec * (i + 1))
    raise RuntimeError(str(last_exc))


def enrichr_to_df(enr, gene_set: str):
    res = enr.results
    if isinstance(res, dict):
        if gene_set in res:
            return res[gene_set].copy()
        for v in res.values():
            if hasattr(v, "copy") and len(v):
                return v.copy()
    return res.copy()


def _col_index(hdr: list[str], *candidates: str) -> int | None:
    """First matching column (case-insensitive), left to right in header."""
    want = [c.strip().lower() for c in candidates]
    for i, h in enumerate(hdr):
        hl = h.strip().lower()
        if hl in want:
            return i
    return None


def plot_top_terms(tsv_path: Path, out_png: Path, title: str, top_n: int, dpi: int) -> None:
    from matplotlib import pyplot as plt

    lines = tsv_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        return
    hdr = lines[0].split("\t")
    term_i = _col_index(hdr, "Term")
    p_i = _col_index(
        hdr,
        "Adjusted P-value",
        "Adjusted P-value (FDR)",
        "Adjusted P-value (Bonferroni)",
        "P-value",
    )
    if term_i is None or p_i is None:
        print(
            "[WARN] TSV missing Term-like or adjusted P-value column; skip plot. "
            f"Columns: {hdr[:12]}{'...' if len(hdr) > 12 else ''}",
            file=sys.stderr,
        )
        return
    top: list[tuple[str, float]] = []
    for line in lines[1:]:
        a = line.split("\t")
        if len(a) <= max(term_i, p_i):
            continue
        term = a[term_i].strip()
        try:
            p = float(a[p_i])
        except (ValueError, IndexError):
            continue
        if term and p > 0 and p <= 1:
            top.append((term, -math.log10(p)))
    top.sort(key=lambda x: x[1], reverse=True)
    top = top[:top_n]
    if not top:
        return
    labels = [t[0][:80] for t in top][::-1]
    vals = [t[1] for t in top][::-1]
    plt.figure(figsize=(10, max(4, 0.35 * len(labels))))
    plt.barh(labels, vals)
    plt.xlabel("-log10(adjusted p-value)")
    plt.title(title)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=dpi)
    plt.close()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", type=Path, required=True, help="Input file")
    p.add_argument(
        "--source",
        choices=("genes", "annotate", "bed"),
        required=True,
        help="genes: one symbol/line; annotate: TSV with gene column; bed: use BED name column as symbol",
    )
    p.add_argument(
        "--gene-column",
        default="Gene Name",
        help='For --source annotate: column name (e.g. "Gene Name") or 1-based column index',
    )
    p.add_argument(
        "--bed-gene-col",
        type=int,
        default=4,
        help="For --source bed: 1-based column index for gene symbol (default 4 = BED name)",
    )
    p.add_argument("--gene-set", default=DEFAULT_GENE_SET, help="Enrichr library name")
    p.add_argument("--organism", default="human", help="gseapy enrichr organism")
    p.add_argument(
        "--out-prefix",
        type=Path,
        required=True,
        help="Writes <prefix>_enrichr.tsv and <prefix>_top_terms.png (unless --no-plot)",
    )
    p.add_argument("--top-n", type=int, default=15, help="Max terms in bar plot")
    p.add_argument("--no-plot", action="store_true", help="TSV only")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--retries", type=int, default=3)
    p.add_argument("--retry-wait-sec", type=int, default=8)
    args = p.parse_args()

    inp = args.input.resolve()
    if not inp.is_file():
        print(f"ERROR: not a file: {inp}", file=sys.stderr)
        return 1

    if args.source == "genes":
        genes = load_genes_plain(inp)
    elif args.source == "annotate":
        genes = load_genes_annotate(inp, args.gene_column)
    else:
        genes = load_genes_bed(inp, args.bed_gene_col)

    print(f"[go_enrichr] unique genes: {len(genes)} from {inp} ({args.source})", file=sys.stderr)
    if not genes:
        print("ERROR: no genes parsed; check --source / column settings.", file=sys.stderr)
        return 1

    prefix = args.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out_tsv = Path(str(prefix) + "_enrichr.tsv")
    out_png = Path(str(prefix) + "_top_terms.png")

    with tempfile.TemporaryDirectory(prefix="goenr_") as tmpd:
        try:
            enr = run_enrichr_with_retry(
                genes,
                args.gene_set,
                args.organism,
                Path(tmpd),
                retries=args.retries,
                retry_wait_sec=args.retry_wait_sec,
            )
            df = enrichr_to_df(enr, args.gene_set)
            df.to_csv(out_tsv, sep="\t", index=False)
        except Exception as e:
            print(f"ERROR: Enrichr failed: {e}", file=sys.stderr)
            return 1

    print(f"Wrote {out_tsv}", file=sys.stderr)
    if not args.no_plot:
        title = f"{args.gene_set}\n({inp.name}, n={len(genes)})"
        plot_top_terms(out_tsv, out_png, title, args.top_n, args.dpi)
        if out_png.is_file():
            print(f"Wrote {out_png}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
