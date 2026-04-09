#!/usr/bin/env python3
"""Enrichr GO BP for Method 2 ABC arms and ER/p65-only control gene lists (MACS + HOMER)."""

from __future__ import annotations

import argparse
import csv
import time
import sys
import tempfile
from pathlib import Path
import math

from matplotlib import pyplot as plt

DEFAULT_GENE_SET = "GO_Biological_Process_2023"


def list_display_name(stem: str) -> str:
    if stem.startswith("genes_"):
        return stem[len("genes_") :]
    return stem


def load_genes_txt(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text().splitlines():
        g = line.strip().upper()
        if g:
            names.append(g)
    return sorted(set(names))


def run_enrichr(genes: list[str], gene_set: str, tmp_parent: Path):
    import gseapy as gp

    tmpd = tmp_parent / "enr"
    tmpd.mkdir(parents=True, exist_ok=True)
    return gp.enrichr(
        gene_list=genes,
        gene_sets=[gene_set],
        organism="human",
        no_plot=True,
        outdir=str(tmpd),
        verbose=False,
    )


def run_enrichr_with_retry(
    genes: list[str],
    gene_set: str,
    tmp_parent: Path,
    retries: int,
    retry_wait_sec: int,
):
    last_exc = None
    for i in range(retries + 1):
        try:
            return run_enrichr(genes, gene_set, tmp_parent)
        except Exception as e:  # network/service throttling from Enrichr
            last_exc = e
            if i >= retries:
                raise
            time.sleep(retry_wait_sec * (i + 1))
    raise RuntimeError(f"unreachable retry flow: {last_exc}")


def enrichr_to_df(enr, gene_set: str):
    res = enr.results
    if isinstance(res, dict):
        if gene_set in res:
            return res[gene_set].copy()
        for v in res.values():
            if hasattr(v, "copy") and len(v):
                return v.copy()
    return res.copy()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cobinding-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent.parent,
        help="analysis/p65_ER_cobinding directory",
    )
    p.add_argument("--gene-set", default=DEFAULT_GENE_SET)
    p.add_argument("--retries", type=int, default=3)
    p.add_argument("--retry-wait-sec", type=int, default=8)
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Default: <cobinding-root>/downstream/go_enrichment",
    )
    args = p.parse_args()
    root = args.cobinding_root.resolve()
    out_dir = (args.out_dir or (root / "downstream/go_enrichment")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = out_dir / "method2_go_manifest.tsv"
    rows_manifest: list[tuple[str, str, int, str]] = []

    with tempfile.TemporaryDirectory(prefix="m2go_") as tmpd:
        tmp_parent = Path(tmpd)
        for caller in ("macs", "homer"):
            ld = root / "downstream/lists" / caller
            if not ld.is_dir():
                print(f"[WARN] missing lists dir {ld}", file=sys.stderr)
                continue
            for txt in sorted(ld.glob("genes_Method2_*.txt")):
                genes = load_genes_txt(txt)
                base = list_display_name(txt.stem)
                out_tsv = out_dir / f"{base}_{caller}_go_bp_enrichr.tsv"
                n_genes = len(genes)
                rows_manifest.append((base, caller, n_genes, str(out_tsv)))
                if n_genes == 0:
                    out_tsv.write_text("")
                    print(f"[skip] 0 genes: {txt.name} ({caller})", file=sys.stderr)
                    continue
                try:
                    enr = run_enrichr_with_retry(
                        genes,
                        args.gene_set,
                        tmp_parent,
                        retries=args.retries,
                        retry_wait_sec=args.retry_wait_sec,
                    )
                    df = enrichr_to_df(enr, args.gene_set)
                    df.to_csv(out_tsv, sep="\t", index=False)
                except Exception as e:
                    print(f"[WARN] Enrichr failed {txt} ({caller}): {e}", file=sys.stderr)
                    out_tsv.write_text("")

    with manifest_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["list_name", "caller", "n_genes", "path"])
        for row in sorted(rows_manifest, key=lambda r: (r[0], r[1])):
            w.writerow(row)

    # Build per-caller, per-arm top-term figures (one PNG per arm).
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    abc_names = ["Method2_A_cross", "Method2_B_cotreat", "Method2_C_cotreat_exclusive"]
    for caller in ("macs", "homer"):
        for abc in abc_names:
            tsv = out_dir / f"{abc}_{caller}_go_bp_enrichr.tsv"
            rows = tsv.read_text().splitlines() if tsv.exists() else []
            if len(rows) < 2:
                continue
            hdr = rows[0].split("\t")
            if "Term" not in hdr or "Adjusted P-value" not in hdr:
                continue
            term_i = hdr.index("Term")
            p_i = hdr.index("Adjusted P-value")
            top: list[tuple[str, float]] = []
            for line in rows[1:]:
                a = line.split("\t")
                if len(a) <= max(term_i, p_i):
                    continue
                term = a[term_i].strip()
                try:
                    p = float(a[p_i])
                except ValueError:
                    continue
                if term and p > 0:
                    top.append((term, -math.log10(p)))
            top.sort(key=lambda x: x[1], reverse=True)
            top = top[:12]
            if not top:
                continue
            labels = [t[0][:70] for t in top][::-1]
            vals = [t[1] for t in top][::-1]
            plt.figure(figsize=(10, 5))
            plt.barh(labels, vals)
            plt.xlabel("-log10(adjusted p-value)")
            plt.title(f"{abc}: GO BP top terms ({caller.upper()})")
            plt.tight_layout()
            plt.savefig(fig_dir / f"{abc}_go_top_terms_{caller}.png", dpi=150)
            plt.close()

    # Remove legacy combined A/B/C figure name if present.
    for caller in ("macs", "homer"):
        (fig_dir / f"method2_abc_go_top_terms_{caller}.png").unlink(missing_ok=True)

    print(f"Wrote {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
