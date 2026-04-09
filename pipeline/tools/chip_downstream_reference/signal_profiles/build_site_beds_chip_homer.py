#!/usr/bin/env python3
"""Build Cotreat_cobind / p65-only / ER-only (sampled) BEDs for ChIP signal profiles (HOMER-based, hg38).

Region sets requested:
  - cobind_E2TNF.bed: Method2 Arm B (cotreat) cobinding bed (HOMER-derived); heatmaps label row "Cotreat_cobind"
  - p65_TNF_only: p65_TNF consensus (HOMER) minus cotreat cobind regions
  - ER_E2_only_10pct: ER_E2 consensus (HOMER) excluding p65_TNF consensus expanded by PROX, then 10% sample
"""

from __future__ import annotations

import argparse
import math
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

_STD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def read_bed(path: Path, require_name: bool = False) -> list[tuple[str, int, int, str]]:
    rows: list[tuple[str, int, int, str]] = []
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0].strip()
            if chrom not in _STD_CHROMS:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            name = parts[3].strip() if len(parts) >= 4 else "."
            if require_name and name == ".":
                continue
            rows.append((chrom, start, end, name))
    return rows


def write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for chrom, s, e, name in rows:
            f.write(f"{chrom}\t{s}\t{e}\t{name}\n")


def sort_bed4_inplace(path: Path) -> None:
    subprocess.run(["sort", "-k1,1", "-k2,2n", "-o", str(path), str(path)], check=True)


def slop_merge_mask(rows: list[tuple[str, int, int, str]], prox: int, genome_sizes: Path, out_mask: Path) -> None:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False, encoding="utf-8") as tmp:
        tmp_path = Path(tmp.name)
        for chrom, s, e, _n in rows:
            tmp.write(f"{chrom}\t{s}\t{e}\n")
    try:
        sort1 = subprocess.run(["sort", "-k1,1", "-k2,2n", str(tmp_path)], check=True, capture_output=True)
        slop = subprocess.run(
            ["bedtools", "slop", "-b", str(prox), "-g", str(genome_sizes), "-i", "stdin"],
            input=sort1.stdout,
            check=True,
            capture_output=True,
        )
        merged = subprocess.run(["bedtools", "merge", "-i", "stdin"], input=slop.stdout, check=True, capture_output=True)
        out_mask.parent.mkdir(parents=True, exist_ok=True)
        out_mask.write_bytes(merged.stdout)
    finally:
        tmp_path.unlink(missing_ok=True)


def subtract_bed(a_rows: list[tuple[str, int, int, str]], mask_bed: Path, out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False, encoding="utf-8") as tmp_a:
        a_tmp = Path(tmp_a.name)
        for chrom, s, e, name in a_rows:
            tmp_a.write(f"{chrom}\t{s}\t{e}\t{name}\n")

    fd, a_sorted_str = tempfile.mkstemp(suffix=".bed", dir=str(out_path.parent))
    os.close(fd)
    a_sorted = Path(a_sorted_str)
    try:
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-o", str(a_sorted), str(a_tmp)], check=True)
        with out_path.open("w") as out_f:
            subprocess.run(["bedtools", "intersect", "-v", "-a", str(a_sorted), "-b", str(mask_bed)], check=True, stdout=out_f)
    finally:
        a_tmp.unlink(missing_ok=True)
        a_sorted.unlink(missing_ok=True)

    return sum(1 for _ in out_path.open())


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cobinding-root", type=Path, default=root, help="analysis/p65_ER_cobinding directory")
    p.add_argument("--genome-sizes", type=Path, default=Path("/mnt/share/archive/bkup/ref/genome/hg38/hg38.chrom.sizes"))
    p.add_argument("--prox", type=int, default=2000)
    p.add_argument("--sample-frac", type=float, default=0.10)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    root = args.cobinding_root.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    cobind_bed_in = root / "downstream/lists/homer/method2_B_cotreat.bed"
    p65_cons_in = root / "homer/method2_A_cross/p65_TNF_consensus.bed"
    er_cons_in = root / "homer/method2_A_cross/ER_E2_consensus.bed"

    for f in (cobind_bed_in, p65_cons_in, er_cons_in, args.genome_sizes):
        if not f.is_file():
            print(f"Missing required input: {f}", file=sys.stderr)
            raise SystemExit(1)

    cobind = read_bed(cobind_bed_in)
    p65_all = read_bed(p65_cons_in, require_name=True)
    er_all = read_bed(er_cons_in, require_name=True)

    cobind_ids = {r[3] for r in cobind if r[3] != "."}
    p65_only = [r for r in p65_all if r[3] not in cobind_ids]

    # ER-only: subtract prox-expanded p65 mask from ER consensus
    mask = out / "_p65_prox_merged.bed"
    slop_merge_mask(p65_all, args.prox, args.genome_sizes, mask)

    er_only_full = out / "_ER_only_full.bed"
    n_er_only = subtract_bed(er_all, mask, er_only_full)

    rng = random.Random(args.seed)
    er_lines = [ln for ln in er_only_full.read_text().splitlines() if ln.strip()]
    if not er_lines:
        chosen = []
    else:
        k = min(len(er_lines), max(1, math.ceil(len(er_lines) * args.sample_frac)))
        chosen = rng.sample(er_lines, k=k)
        chosen.sort(key=lambda ln: (ln.split("\t")[0], int(ln.split("\t")[1])))

    cobind_out = out / "cobind_E2TNF.bed"
    p65_only_out = out / "p65_TNF_only.bed"
    er_only_out = out / "ER_E2_only_10pct.bed"

    write_bed(cobind_out, cobind)
    write_bed(p65_only_out, p65_only)
    er_only_out.write_text("\n".join(chosen) + ("\n" if chosen else ""))

    qc = out / "site_bed_qc.txt"
    qc.write_text(
        "\n".join(
            [
                f"cobind_bed\t{cobind_bed_in}",
                f"p65_consensus\t{p65_cons_in}",
                f"er_consensus\t{er_cons_in}",
                f"genome_sizes\t{args.genome_sizes}",
                f"prox_bp\t{args.prox}",
                f"sample_frac\t{args.sample_frac}",
                f"seed\t{args.seed}",
                f"n_cobind\t{len(cobind)}",
                f"n_p65_consensus_stdChrom\t{len(p65_all)}",
                f"n_p65_only\t{len(p65_only)}",
                f"n_er_consensus_stdChrom\t{len(er_all)}",
                f"n_er_only_before_sample\t{n_er_only}",
                f"n_er_only_sampled\t{len(chosen)}",
                "",
            ]
        )
    )

    mask.unlink(missing_ok=True)
    er_only_full.unlink(missing_ok=True)
    print(qc.read_text())


if __name__ == "__main__":
    main()

