#!/usr/bin/env python3
"""
dupescan.py -- systematic cross-volume duplicate-file investigation.

Four phases, run in order. Each writes into OUTDIR and is independently
resumable, so a scan interrupted by an NFS stall picks up where it stopped.

  index    metadata-only walk of every enabled root in roots.tsv.
           No file contents are read.  ->  index/<label>.tsv.gz
  pair     pure-CPU grouping of the indexes by size, after collapsing
           hardlinks and alias mounts on (dev,ino).  ->  candidates.tsv.gz
  verify   reads only head+tail 64 KiB of each candidate and fingerprints it,
           which is what separates real duplicates from same-size coincidences.
           --full hashes candidates end to end instead.  ->  fp_cache.tsv
  report   cross-root byte matrix, largest duplicate groups, and the
           directory-containment rollup.  ->  report/

The directory rollup is the phase that answers the actual question. Individual
duplicate pairs number in the millions and are not actionable; "this tree is
97% contained in that tree" is.

Nothing in this script deletes, moves, or writes to any scanned root.
"""

import argparse
import errno
import gzip
import hashlib
import os
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOTS = os.path.join(HERE, "roots.tsv")
DEFAULT_OUT = os.path.join(HERE, "out")

# Synology / NFS bookkeeping trees. Indexing them produces noise, not signal.
PRUNE_NAMES = ["@eaDir", "#recycle", "#snapshot", ".snapshot", ".zfs",
               "@tmp", "lost+found"]

EDGE = 65536          # bytes hashed from each end during verify
MIN_SIZE = 4096       # files below this are ignored; they cannot add up
WORKERS = 8           # concurrent readers; NAS is the bottleneck, not CPU


def log(msg):
    sys.stderr.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))
    sys.stderr.flush()


def human(n):
    n = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
        if abs(n) < 1024.0:
            return "%.1f %s" % (n, unit)
        n /= 1024.0
    return "%.1f EiB" % n


def read_roots(path):
    roots = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            label, root, enabled = parts[0], parts[1], parts[2]
            roots.append((label, root, enabled == "1"))
    return roots


# --------------------------------------------------------------------------
# phase: index
# --------------------------------------------------------------------------

def mount_of(path):
    """(mount_point, device_spec) of the filesystem holding `path`."""
    path = os.path.realpath(path)
    best = ("", "")
    try:
        with open("/proc/mounts") as fh:
            for line in fh:
                f = line.split()
                if len(f) < 2:
                    continue
                dev, mp = f[0], f[1].replace("\\040", " ")
                if (path == mp or path.startswith(mp.rstrip("/") + "/")) \
                        and len(mp) > len(best[0]):
                    best = (mp, dev)
    except OSError:
        pass
    return best


def check_alias_mounts(roots):
    """Refuse to scan two paths that are the same export mounted twice.

    linux01 has several of these (the same NFS export mounted at two paths,
    sometimes via different server IPs). Indexing both makes every file its
    own duplicate and silently doubles every byte total, so this is a hard
    error rather than a warning.
    """
    by_dev = defaultdict(list)
    for label, root, enabled in roots:
        if not enabled:
            continue
        mp, dev = mount_of(root)
        if dev:
            by_dev[dev].append((label, root, mp))
    bad = []

    # Second, independent signal: two different servers can export the same
    # backing filesystem under different hostnames, which /proc/mounts cannot
    # see. Identical root-directory inode + size gives it away.
    by_stat = defaultdict(list)
    for label, root, enabled in roots:
        if not enabled:
            continue
        try:
            st = os.stat(root)
        except OSError:
            continue
        by_stat[(st.st_ino, st.st_size, st.st_mtime)].append((label, root))
    for _k, entries in by_stat.items():
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                bad.append(("identical root inode", entries[i][0], entries[i][1],
                            entries[j][0], entries[j][1]))

    for dev, entries in by_dev.items():
        if len(entries) < 2:
            continue
        # Same export is only a real alias if the scanned subtrees overlap.
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                la, ra, ma = entries[i]
                lb, rb, mb = entries[j]
                sa = os.path.realpath(ra)[len(ma):] or "/"
                sb = os.path.realpath(rb)[len(mb):] or "/"
                if sa.startswith(sb) or sb.startswith(sa):
                    bad.append((dev, la, ra, lb, rb))
    bad = sorted({(d, la, ra, lb, rb) for d, la, ra, lb, rb in bad},
                 key=lambda t: (t[1], t[3]))
    if bad:
        sys.stderr.write("\nERROR: roots.tsv enables the same export twice.\n")
        sys.stderr.write("Every file below would be counted as its own duplicate.\n\n")
        for dev, la, ra, lb, rb in bad:
            sys.stderr.write("  %s\n    %-14s %s\n    %-14s %s\n"
                             % (dev, la, ra, lb, rb))
        sys.stderr.write("\nDisable one of each pair (set its enabled column to 0).\n")
        sys.exit(2)


def phase_index(args):
    roots = read_roots(args.roots)
    check_alias_mounts(roots)
    idx_dir = os.path.join(args.outdir, "index")
    os.makedirs(idx_dir, exist_ok=True)

    for label, root, enabled in roots:
        if not enabled:
            continue
        if args.only and label not in args.only:
            continue
        dest = os.path.join(idx_dir, label + ".tsv.gz")
        if os.path.exists(dest) and not args.force:
            log("index %-14s SKIP (exists; --force to redo)" % label)
            continue
        if not os.path.isdir(root):
            log("index %-14s SKIP (not a directory: %s)" % (label, root))
            continue

        # Probe first: a dead volume must fail here, not 3 hours in.
        try:
            os.listdir(root)
        except OSError as exc:
            log("index %-14s FAIL (%s) -- root unreadable, skipping"
                % (label, errno.errorcode.get(exc.errno, exc.errno)))
            continue

        cmd = ["find", root, "-xdev"]
        cmd += ["("] + sum([["-name", n, "-o"] for n in PRUNE_NAMES], [])[:-1] \
             + [")", "-prune", "-o"]
        cmd += ["-type", "f", "-printf", r"%s\t%D\t%i\t%n\t%T@\t%p\0"]
        if _has("ionice"):
            cmd = ["ionice", "-c3", "-t"] + cmd
        cmd = ["nice", "-n", "19"] + cmd

        log("index %-14s scanning %s" % (label, root))
        t0 = time.time()
        tmp = dest + ".part"
        n = 0
        errs = 0
        with gzip.open(tmp, "wt", encoding="utf-8", errors="surrogateescape") as out:
            out.write("size\tdev\tino\tnlink\tmtime\tpath\n")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            buf = b""
            while True:
                chunk = proc.stdout.read(1 << 20)
                if not chunk:
                    break
                buf += chunk
                recs = buf.split(b"\0")
                buf = recs.pop()
                for rec in recs:
                    if not rec:
                        continue
                    out.write(rec.decode("utf-8", "surrogateescape") + "\n")
                    n += 1
            stderr = proc.stderr.read()
            proc.wait()
            errs = stderr.count(b"\n")
        os.rename(tmp, dest)
        dt = time.time() - t0
        log("index %-14s %d files in %.0fs (%.0f files/s, %d fs errors)"
            % (label, n, dt, n / dt if dt else 0, errs))


def _has(prog):
    for d in os.environ.get("PATH", "").split(os.pathsep):
        if os.access(os.path.join(d, prog), os.X_OK):
            return True
    return False


def load_index(idx_dir, labels=None):
    """Yield (label, size, dev, ino, nlink, mtime, path)."""
    for fn in sorted(os.listdir(idx_dir)):
        if not fn.endswith(".tsv.gz"):
            continue
        label = fn[:-7]
        if labels and label not in labels:
            continue
        with gzip.open(os.path.join(idx_dir, fn), "rt",
                       encoding="utf-8", errors="surrogateescape") as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t", 5)
                if len(f) < 6:
                    continue
                yield (label, int(f[0]), int(f[1]), int(f[2]),
                       int(f[3]), f[4], f[5])


# --------------------------------------------------------------------------
# phase: pair
# --------------------------------------------------------------------------

def phase_pair(args):
    idx_dir = os.path.join(args.outdir, "index")
    by_size = defaultdict(list)
    seen_inode = {}
    stats = defaultdict(int)

    for label, size, dev, ino, nlink, mtime, path in load_index(idx_dir):
        stats["files"] += 1
        stats["bytes"] += size
        if size < args.min_size:
            stats["skipped_small"] += 1
            continue
        key = (dev, ino)
        if key in seen_inode:
            # Same physical file reached twice: a hardlink, or an alias mount
            # that slipped into roots.tsv. Either way it is not a duplicate.
            stats["collapsed_inode"] += 1
            continue
        seen_inode[key] = 1
        by_size[size].append((label, dev, ino, path))

    dest = os.path.join(args.outdir, "candidates.tsv.gz")
    ngroups = ncand = cand_bytes = 0
    with gzip.open(dest, "wt", encoding="utf-8", errors="surrogateescape") as out:
        out.write("size\tlabel\tdev\tino\tpath\n")
        for size, members in by_size.items():
            if len(members) < 2:
                continue
            ngroups += 1
            for label, dev, ino, path in members:
                out.write("%d\t%s\t%d\t%d\t%s\n" % (size, label, dev, ino, path))
                ncand += 1
            cand_bytes += size * (len(members) - 1)

    log("pair  indexed %d files / %s" % (stats["files"], human(stats["bytes"])))
    log("pair  skipped %d files < %s" % (stats["skipped_small"], human(args.min_size)))
    log("pair  collapsed %d hardlink/alias-mount inodes" % stats["collapsed_inode"])
    log("pair  %d same-size groups, %d candidate files" % (ngroups, ncand))
    log("pair  UPPER BOUND on reclaimable space: %s (before content check)"
        % human(cand_bytes))
    log("pair  -> %s" % dest)


# --------------------------------------------------------------------------
# phase: verify
# --------------------------------------------------------------------------

def fingerprint(path, size, full=False):
    h = hashlib.blake2b(digest_size=16)
    h.update(b"%d|" % size)
    try:
        with open(path, "rb", buffering=0) as fh:
            if full:
                while True:
                    b = fh.read(1 << 22)
                    if not b:
                        break
                    h.update(b)
            else:
                h.update(fh.read(EDGE))
                if size > 2 * EDGE:
                    fh.seek(-EDGE, os.SEEK_END)
                    h.update(fh.read(EDGE))
    except OSError as exc:
        return "ERR:" + errno.errorcode.get(exc.errno, str(exc.errno))
    return h.hexdigest()


def phase_verify(args):
    cand = os.path.join(args.outdir, "candidates.tsv.gz")
    cache_path = os.path.join(args.outdir,
                              "fp_full.tsv" if args.full else "fp_cache.tsv")

    done = {}
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8", errors="surrogateescape") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t", 1)
                if len(f) == 2:
                    done[f[0]] = f[1]
        log("verify resuming with %d cached fingerprints" % len(done))

    todo = []
    deferred = [0]
    with gzip.open(cand, "rt", encoding="utf-8", errors="surrogateescape") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t", 4)
            if len(f) < 5:
                continue
            size, dev, ino, path = int(f[0]), f[2], f[3], f[4]
            if size < args.min_size:
                deferred[0] += 1
                continue
            key = "%s:%s" % (dev, ino)
            if key not in done:
                todo.append((key, path, size))

    log("verify %d files to fingerprint (%s mode)"
        % (len(todo), "FULL" if args.full else "head+tail"))
    if deferred[0]:
        log("verify %d candidates deferred below --min-size %s; rerun with a "
            "lower floor to cover them" % (deferred[0], human(args.min_size)))
    if not todo:
        log("verify nothing to do")
        return

    t0 = time.time()
    n = 0
    with open(cache_path, "a", encoding="utf-8", errors="surrogateescape") as out:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {}
            it = iter(todo)
            # Bounded in-flight window: never materialise millions of futures.
            def submit_next():
                try:
                    key, path, size = next(it)
                except StopIteration:
                    return False
                futs[pool.submit(fingerprint, path, size, args.full)] = key
                return True

            for _ in range(args.workers * 8):
                if not submit_next():
                    break
            while futs:
                from concurrent.futures import wait, FIRST_COMPLETED
                ready, _pending = wait(list(futs), return_when=FIRST_COMPLETED)
                for fut in ready:
                    key = futs.pop(fut)
                    out.write("%s\t%s\n" % (key, fut.result()))
                    n += 1
                    if n % 5000 == 0:
                        out.flush()
                        rate = n / (time.time() - t0)
                        log("verify %d/%d (%.0f files/s, eta %.0f min)"
                            % (n, len(todo), rate,
                               (len(todo) - n) / rate / 60 if rate else 0))
                    submit_next()
    log("verify done: %d fingerprints in %.0fs -> %s"
        % (n, time.time() - t0, cache_path))


# --------------------------------------------------------------------------
# phase: report
# --------------------------------------------------------------------------

def trunc_dir(path, root, depth):
    rel = os.path.relpath(os.path.dirname(path), root)
    if rel == ".":
        return root
    parts = rel.split(os.sep)[:depth]
    return os.path.join(root, *parts)


def phase_report(args):
    roots = {l: r for l, r, _ in read_roots(args.roots)}
    cache_path = os.path.join(args.outdir,
                              "fp_full.tsv" if args.full else "fp_cache.tsv")
    fp = {}
    with open(cache_path, encoding="utf-8", errors="surrogateescape") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t", 1)
            if len(f) == 2:
                fp[f[0]] = f[1]

    groups = defaultdict(list)
    nunverified = 0
    unverified_bytes = [0]
    cand = os.path.join(args.outdir, "candidates.tsv.gz")
    nerr = 0
    with gzip.open(cand, "rt", encoding="utf-8", errors="surrogateescape") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t", 4)
            if len(f) < 5:
                continue
            size, label, dev, ino, path = int(f[0]), f[1], f[2], f[3], f[4]
            d = fp.get("%s:%s" % (dev, ino))
            if d is None:
                nunverified += 1
                unverified_bytes[0] += size
                continue
            if d.startswith("ERR:"):
                nerr += 1
                continue
            groups[(size, d)].append((label, path))

    rep = os.path.join(args.outdir, "report")
    os.makedirs(rep, exist_ok=True)

    matrix = defaultdict(lambda: [0, 0])        # (rootA,rootB) -> [files, bytes]
    within = defaultdict(lambda: [0, 0])        # rootA -> [files, bytes]
    dirpair = defaultdict(lambda: [0, 0])       # (lA,dA,lB,dB) -> [files, bytes]
    top = []
    capped = [0]
    total_waste = 0
    ngroups = 0

    for (size, _d), members in groups.items():
        if len(members) < 2:
            continue
        ngroups += 1
        waste = size * (len(members) - 1)
        total_waste += waste
        top.append((waste, size, len(members), members[0][1], members[1][1]))

        labels = sorted({m[0] for m in members})
        if len(labels) == 1:
            within[labels[0]][0] += len(members) - 1
            within[labels[0]][1] += waste
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                cell = matrix[(labels[i], labels[j])]
                cell[0] += 1
                cell[1] += size

        # One representative per rolled-up directory, so 40 copies sitting in
        # the same directory do not swamp the tree-vs-tree signal. Pairing is
        # by directory rather than by root, so duplication *inside* one root
        # (the usual shape of a half-finished restore) shows up too.
        reps = {}
        for label, path in members:
            d = trunc_dir(path, roots.get(label, ""), args.depth)
            reps.setdefault((label, d), path)
        rl = sorted(reps)
        if len(rl) > args.max_dirs_per_group:
            capped[0] += 1
            rl = rl[:args.max_dirs_per_group]
        for i in range(len(rl)):
            for j in range(i + 1, len(rl)):
                (la, da), (lb, db) = rl[i], rl[j]
                cell = dirpair[(la, da, lb, db)]
                cell[0] += 1
                cell[1] += size

    # Directory totals, for containment percentages.
    idx_dir = os.path.join(args.outdir, "index")
    dtot = defaultdict(lambda: [0, 0])
    for label, size, dev, ino, nlink, mtime, path in load_index(idx_dir):
        d = trunc_dir(path, roots.get(label, ""), args.depth)
        dtot[(label, d)][0] += 1
        dtot[(label, d)][1] += size

    with open(os.path.join(rep, "root_matrix.tsv"), "w") as out:
        out.write("root_a\troot_b\tdup_files\tdup_bytes\tdup_human\n")
        for (a, b), (n, by) in sorted(matrix.items(), key=lambda kv: -kv[1][1]):
            out.write("%s\t%s\t%d\t%d\t%s\n" % (a, b, n, by, human(by)))
        for a, (n, by) in sorted(within.items(), key=lambda kv: -kv[1][1]):
            out.write("%s\t(same root)\t%d\t%d\t%s\n" % (a, n, by, human(by)))

    top.sort(reverse=True)
    with open(os.path.join(rep, "top_groups.tsv"), "w",
              encoding="utf-8", errors="surrogateescape") as out:
        out.write("wasted_bytes\twasted_human\tfile_size\tn_copies\tpath_a\tpath_b\n")
        for waste, size, n, pa, pb in top[:args.top]:
            out.write("%d\t%s\t%d\t%d\t%s\t%s\n" % (waste, human(waste), size, n, pa, pb))

    rows = []
    for (la, da, lb, db), (n, by) in dirpair.items():
        ta = dtot.get((la, da), [0, 0])
        tb = dtot.get((lb, db), [0, 0])
        rows.append((by, n, la, da, ta[1], 100.0 * by / ta[1] if ta[1] else 0.0,
                     lb, db, tb[1], 100.0 * by / tb[1] if tb[1] else 0.0))
    rows.sort(reverse=True)
    with open(os.path.join(rep, "dir_pairs.tsv"), "w",
              encoding="utf-8", errors="surrogateescape") as out:
        out.write("dup_bytes\tdup_human\tn_files\troot_a\tdir_a\tdir_a_bytes\tpct_of_a"
                  "\troot_b\tdir_b\tdir_b_bytes\tpct_of_b\n")
        for by, n, la, da, ab, pa, lb, db, bb, pb in rows[:args.top]:
            out.write("%d\t%s\t%d\t%s\t%s\t%d\t%.1f\t%s\t%s\t%d\t%.1f\n"
                      % (by, human(by), n, la, da, ab, pa, lb, db, bb, pb))

    lines = []
    lines.append("dupescan report  (%s fingerprint)"
                 % ("FULL-content" if args.full else "head+tail+size"))
    lines.append("=" * 62)
    lines.append("duplicate groups        : %d" % ngroups)
    lines.append("unreadable candidates   : %d" % nerr)
    if nunverified:
        lines.append("UNVERIFIED candidates   : %d (%s of candidate bytes) -- "
                     "not counted below; run verify with a lower --min-size"
                     % (nunverified, human(unverified_bytes[0])))
    if capped[0]:
        lines.append("groups capped at %d dirs : %d  (dir rollup only; byte "
                     "totals above are complete)" % (args.max_dirs_per_group, capped[0]))
    lines.append("reclaimable if collapsed: %s" % human(total_waste))
    lines.append("")
    lines.append("By root pair (largest first):")
    for (a, b), (n, by) in sorted(matrix.items(), key=lambda kv: -kv[1][1])[:20]:
        lines.append("  %-14s <-> %-14s %10s  (%d groups)" % (a, b, human(by), n))
    lines.append("")
    lines.append("Within a single root (internal duplication):")
    for a, (n, by) in sorted(within.items(), key=lambda kv: -kv[1][1])[:20]:
        lines.append("  %-14s %10s  (%d redundant files)" % (a, human(by), n))
    lines.append("")
    lines.append("Largest duplicated directory pairs:")
    for by, n, la, da, ab, pa, lb, db, bb, pb in rows[:20]:
        lines.append("  %9s  %s" % (human(by), da))
        lines.append("             = %.0f%% of it, also at %s (%.0f%% of that)"
                     % (pa, db, pb))
    text = "\n".join(lines)
    with open(os.path.join(rep, "summary.txt"), "w",
              encoding="utf-8", errors="surrogateescape") as out:
        out.write(text + "\n")
    print(text)
    log("report -> %s" % rep)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("phase", choices=["index", "pair", "verify", "report"])
    p.add_argument("--roots", default=DEFAULT_ROOTS)
    p.add_argument("--outdir", default=DEFAULT_OUT)
    p.add_argument("--only", nargs="*", help="restrict index phase to these labels")
    p.add_argument("--force", action="store_true", help="re-index roots already done")
    p.add_argument("--min-size", type=int, default=MIN_SIZE)
    p.add_argument("--workers", type=int, default=WORKERS)
    p.add_argument("--full", action="store_true",
                   help="verify/report on whole-file hashes instead of head+tail")
    p.add_argument("--depth", type=int, default=4,
                   help="directory rollup depth below each root (report)")
    p.add_argument("--top", type=int, default=200)
    p.add_argument("--max-dirs-per-group", type=int, default=12,
                   help="cap on distinct directories paired per duplicate group")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    {"index": phase_index, "pair": phase_pair,
     "verify": phase_verify, "report": phase_report}[args.phase](args)


if __name__ == "__main__":
    main()
