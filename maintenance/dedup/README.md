# scripts/maintenance/dedup — cross-volume duplicate investigation

Read-only survey of which files exist in more than one place across the
mounted volumes. Nothing here deletes, moves, or writes to a scanned root;
the output is evidence for a human decision.

Motivation: after a drive failure and partial restores, the same data tends to
survive in several places at once — a half-finished copy, a "junk"/"old" tree
someone kept "just in case", an offsite mirror of a volume that has since been
rebuilt. Finding that costs nothing but a metadata walk plus a few kilobytes
read per candidate file.

## Method

Four phases, each resumable, each cheaper than the naive approach:

| Phase | Reads | Cost | Output |
|---|---|---|---|
| `index` | metadata only (`find -printf`) | ~2k files/s per cold NFS export | `out/index/<label>.tsv.gz` |
| `pair` | nothing (pure CPU on the indexes) | seconds | `out/candidates.tsv.gz` |
| `verify` | head + tail 64 KiB of candidates only | ~250 files/s | `out/fp_cache.tsv` |
| `report` | nothing | seconds | `out/report/` |

The point of the split is that **no whole file is ever read**. Two files can
only be identical if they have the same size, so `pair` reduces millions of
files to the few that share a size with something; `verify` then reads 128 KiB
of each. For sequencing data — where same-size files are usually genuinely the
same file — `size + first 64 KiB + last 64 KiB` is decisive in practice. Pass
`--full` to `verify` and `report` when you want whole-file hashes instead;
it is the same pipeline, just orders of magnitude more I/O.

### What the phases guard against

- **Alias mounts.** The same NFS export mounted at two paths would make every
  file its own duplicate and silently double every byte total. `index` refuses
  to start if two enabled roots resolve to one export, detected both from
  `/proc/mounts` and from identical root-directory inodes (which catches the
  case of one filesystem served under two hostnames).
- **Hardlinks.** Collapsed on `(dev, ino)` in `pair`. A hardlink occupies one
  copy of the data, so counting it as reclaimable would overstate the win.
- **Small files.** Ignored below `--min-size` (default 4 KiB). They are the
  overwhelming majority of files and a negligible fraction of the bytes.
- **Bookkeeping trees.** `@eaDir`, `#recycle`, `.snapshot` and friends are
  pruned during the walk.

## Running it

```bash
cd ~/work/scripts/maintenance/dedup
cp roots.tsv.example roots.tsv     # then edit: which volumes, which enabled
PY=~/miniforge3/envs/bio/bin/python3

$PY dupescan.py index            # hours; nice/ionice'd; safe to interrupt
$PY dupescan.py pair             # prints the upper bound on reclaimable space
$PY dupescan.py verify           # confirms which candidates are real
$PY dupescan.py report --depth 4
```

`index` skips roots it has already done, so an interrupted scan resumes by
re-running the same command (`--force` to redo one). `verify` caches every
fingerprint by `dev:ino` and resumes the same way. Scan a subset with
`--only <label> <label>`.

Run `index` one or two roots at a time on a busy NAS. It is `nice -n 19` and
`ionice -c3`, but it is still a few thousand metadata ops per second, and this
NAS serves everyone's home directory.

## Reading the output

`report/summary.txt` is the whole story; the TSVs beside it are the same data
unabridged.

- **By root pair** — bytes duplicated between each pair of volumes. This is
  the "do I have two copies of this volume" number.
- **Within a single root** — redundancy that has accumulated inside one
  volume. On a system that has been through restores, this is usually the
  larger number.
- **Largest duplicated directory pairs** (`dir_pairs.tsv`) — the actionable
  one. Individual duplicate pairs number in the millions and no one can act on
  them; `dir A is 89% contained in dir B` is a decision someone can make in a
  few seconds. `pct_of_a` and `pct_of_b` are each side's duplicated bytes as a
  fraction of that directory's total, so `100% / 40%` reads as "A is entirely
  inside B, and is 40% of B".

A directory pair at 100% is a copy. One at 50–90% is usually a version bump
(`star/2.7.1a` vs `star/2.7.3a`) where the shared files are the inputs and the
differing ones are the built index — real duplication, but deleting either
side loses something.

### The headline number is not the reclaimable number

`reclaimable if collapsed` counts every duplicate group, and **a backup is a
duplicate by construction**. Duplication between a live volume and its offsite
mirror is the backup working correctly; collapsing it would destroy the very
redundancy it exists to provide. Read the two blocks separately:

- *By root pair* involving the offsite mirror — expected, leave alone. What is
  worth checking here is *what* is being mirrored: paying offsite capacity to
  preserve a trash directory is a real finding, and deleting on the live side
  is what fixes it.
- *Within a single root* — one volume holding the same bytes twice. This is
  the unambiguous waste, and usually the smaller of the two numbers.

Anything reported here is a claim about bytes, not about references. A tree
being byte-identical to another says nothing about whether a script, symlink,
or someone's notebook still points at that path.

## Deleting anything

Out of scope for this tool, deliberately. `report` tells you which trees to
look at; confirming a tree is safe to remove needs `--full` verification and a
human who knows what the tree is for. Two files being byte-identical today
says nothing about whether something still points at both paths.
