#!/usr/bin/env python3
"""frag_hole_check.py — per-library check for the read-length hole in paired-end fragment lengths (bio env: pysam, numpy, pandas).

Why: when adapter trimming needs more overlap than a fragment's adapter tail (Trim Galore --stringency 5), fragments 1–4 bp
shorter than the read keep 1–4 bp of adapter on both reads. The mates then overshoot each other: without bowtie2 --dovetail the
pair is not proper and drops out (a hole at read length − 4 … − 1), with it the pair is recorded 2–8 bp too long (a spike at
read length … + 4). scripts/pipeline/cnr trims at stringency 1 since 2026-11-01, which closes it; this check catches a library
where it does not (other adapters, other assays, reads under 50 bp that stay on 5). Read a dip or spike at the read length as
this artefact, not biology.

Fragments are counted as in seq/cnr/CnR_261006_ERa-ARNT-p65_MCF7-HEK293T_Steven/analysis/fragment_size_profiles_v1 (01_frag_hist):
proper pairs, read 1 only, no unmapped / secondary / supplementary, both mates on one chromosome, chromosomes 1–N, X and Y,
deduplicated by fragment start + length, fragments overlapping the blacklist left out, no MAPQ filter. L = the longest read 1
seen. With d = fragments per 10,000 at each length:
  hole  = mean d at L − 4 … L − 1  /  mean d at the flanks L − 9 … L − 5 and L + 5 … L + 9   (about 1 = no hole)
  spike = mean d at L … L + 4      /  the same flanks                                     (about 1 = no spike)
Status: ok, HOLE (hole < --hole-min), SPIKE (spike > --spike-max), both, or LOW_COUNT when the flanks hold fewer than
--min-flank fragments in all. It warns and exits 0; it never stops a pipeline.
The ratios assume the length distribution is smooth across L ± 9. When a real mode sits there, the spike ratio moves without
any artefact: at PE151 the mono-nucleosome peak of histone marks lies at the read length (GS37 H3K27me3, stringency 5 and no
--dovetail, so no remnant pairs kept: spike 1.22–1.25, hole 0.02). Read a histone library's flag against its histogram.

Output: one TSV row per BAM, with a MultiQC custom-content header (name it *_mqc.tsv to have MultiQC show it).

Usage:
  frag_hole_check.py --out align/bam/bamqc/frag_hole_check_mqc.tsv --blacklist hg38-blacklist.v2.bed --jobs 8 align/bam/*_sorted.bam
"""
from __future__ import annotations

import argparse
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

HMAX = 1000                                     # lengths 1..HMAX one by one; longer ones are not needed here
MAIN = re.compile(r"^(chr)?([0-9]+|X|Y)$")
BLACKLIST = None                                # chrom -> (starts, ends), set in main(), inherited by the workers


def read_blacklist(path):
    if not path:
        return {}
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"], comment="#")
    return {c: (g.s.to_numpy(np.int64), g.e.to_numpy(np.int64)) for c, g in b.sort_values(["c", "s"]).groupby("c")}


def count(bam_path):
    """length histogram of deduplicated proper-pair fragments, and the longest read 1."""
    h = np.zeros(HMAX + 2, dtype=np.int64)       # index = length; HMAX + 1 = longer
    read_len = 0
    with pysam.AlignmentFile(bam_path) as bam:
        for chrom in bam.references:
            if not MAIN.match(chrom):
                continue
            st, ln = [], []
            for r in bam.fetch(chrom):
                f = r.flag
                if (f & 0x42) != 0x42 or (f & 0x904) or r.next_reference_id != r.reference_id:
                    continue
                t = r.template_length
                if t == 0:
                    continue
                if r.query_length > read_len:
                    read_len = r.query_length
                st.append(r.reference_start if t > 0 else r.next_reference_start)
                ln.append(abs(t))
            if not st:
                continue
            st, ln = np.asarray(st, np.int64), np.asarray(ln, np.int64)
            key = np.unique(st * 1_000_000 + np.minimum(ln, 999_999))
            st, ln = key // 1_000_000, key % 1_000_000
            if chrom in BLACKLIST:
                bs, be = BLACKLIST[chrom]
                i = np.searchsorted(bs, st + ln, side="left") - 1
                hit = (i >= 0) & (be[np.clip(i, 0, None)] > st)
                st, ln = st[~hit], ln[~hit]
            h += np.bincount(np.minimum(ln, HMAX + 1), minlength=HMAX + 2)
    return bam_path, h, read_len


def score(h, L, hole_min, spike_max, min_flank):
    c = h[1:HMAX + 1].astype(float)              # c[i - 1] = fragments of length i
    total = h.sum()
    row = dict(read_len=L, fragments=int(total), flank_fragments=0, hole=np.nan, spike=np.nan)
    if L < 10 or L + 9 > HMAX or total == 0:
        row["status"] = "NO_DATA"
        return row
    d = 1e4 * c / total
    at = lambda a, b: d[a - 1:b]                 # lengths a..b inclusive
    flank = np.r_[at(L - 9, L - 5), at(L + 5, L + 9)].mean()
    row["flank_fragments"] = int(np.r_[c[L - 10:L - 5], c[L + 4:L + 9]].sum())
    if flank > 0:
        row["hole"], row["spike"] = at(L - 4, L - 1).mean() / flank, at(L, L + 4).mean() / flank
    if row["flank_fragments"] < min_flank or flank == 0:
        row["status"] = "LOW_COUNT"
    else:
        bad = [s for s, v in (("HOLE", row["hole"] < hole_min), ("SPIKE", row["spike"] > spike_max)) if v]
        row["status"] = "+".join(bad) if bad else "ok"
    return row


def main():
    global BLACKLIST
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bams", nargs="+", help="coordinate-sorted, indexed paired-end BAMs")
    ap.add_argument("--out", required=True, help="output TSV (name it *_mqc.tsv for MultiQC)")
    ap.add_argument("--blacklist", default="", help="BED of regions to leave out (optional)")
    ap.add_argument("--jobs", type=int, default=4, help="BAMs counted in parallel (one CPU each)")
    ap.add_argument("--hole-min", type=float, default=0.8)
    ap.add_argument("--spike-max", type=float, default=1.2)
    ap.add_argument("--min-flank", type=int, default=1000, help="fewest fragments in the 10 flank lengths to score")
    ap.add_argument("--suffix", default="_sorted.bam", help="stripped from the BAM name to give the library name")
    a = ap.parse_args()
    for b in a.bams:
        if not (Path(b).is_file() and (Path(b + ".bai").is_file() or Path(b[:-4] + ".bai").is_file())):
            sys.exit(f"[frag_hole_check] missing BAM or index: {b}")
    BLACKLIST = read_blacklist(a.blacklist)
    with Pool(max(1, min(a.jobs, len(a.bams)))) as pool:
        res = pool.map(count, a.bams)
    rows = []
    for b, h, L in res:
        name = Path(b).name
        name = name[:-len(a.suffix)] if a.suffix and name.endswith(a.suffix) else Path(b).stem
        rows.append(dict(library=name, **score(h, L, a.hole_min, a.spike_max, a.min_flank)))
    T = pd.DataFrame(rows)[["library", "read_len", "fragments", "flank_fragments", "hole", "spike", "status"]]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    head = (
        "# id: 'frag_hole_check'\n"
        "# section_name: 'Read-length hole'\n"
        f"# description: 'Fragment-length density at read length − 4 … − 1 (hole) and read length … + 4 (spike), each over the flanks; "
        f"about 1 = no artefact. Flagged when hole < {a.hole_min} or spike > {a.spike_max}: adapter remnants on fragments just shorter "
        f"than the read (scripts/pipeline/tools/frag_hole_check.py).'\n"
        "# plot_type: 'table'\n"
    )
    out.write_text(head + T.to_csv(sep="\t", index=False, float_format="%.4f"))
    for r in T.itertuples():
        msg = f"{r.library}: read {r.read_len} bp, hole {r.hole:.3f}, spike {r.spike:.3f} ({r.status})"
        print(("[frag_hole_check] " if r.status == "ok" else "[WARN] frag_hole_check ") + msg, flush=True)
    print(f"[frag_hole_check] {len(T)} libraries, {int((T.status != 'ok').sum())} flagged -> {out}")


if __name__ == "__main__":
    main()
