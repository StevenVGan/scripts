# `scripts/env/` — conda environment tracking

This directory tracks the conda environments the workspace depends on, so old
projects stay reproducible after upgrades. `bio` is the workhorse for the bulk
pipelines (cutrun, csRNA, proseq, atac) and the shared tools; `sc`, `meth`,
`rna`, `primer`, `pwm2` and `pwm2-pb` serve single-cell, methylome, RNA-seq,
primer-design and PWM work respectively.

## Files — three kinds, updated differently

Telling them apart matters: **a `conda env export` over a spec destroys it.**
The marker is pinning: an export pins every dependency `name=version=build` and
usually ends with a `prefix:` line; a spec pins almost nothing. (Don't rely on
`prefix:` alone — `pwm2.yml` is an export that had its `prefix:` line stripped.)

- **Spec** — `sc.yml`, `meth.yml`. Hand-written, grouped, commented, loosely
  pinned (usually only `python=`). Records *intent* and must stay solvable on a
  fresh box. **Hand-edit these; never `conda env export` over them.** Add the
  package plus a one-line comment saying what needs it. `meth.yml`, for example,
  carries the reasoning for pinning `setuptools<81` and for replacing
  `ucsc-liftover` with `CrossMap` on glibc 2.23 — an export would erase it.
- **Export** — `bio.yml`, `rna.yml`, `primer.yml`, `pwm2.yml`, `pwm2-pb.yml`.
  `conda env export` output: every dependency pinned `name=version=build`, plus
  a `prefix:` line. Re-export after any install.
- **Lock** — `lock/<env>.YYYY-MM-DD.yml`. Always a full export, for *both* kinds
  above. The lock is where exactness lives, which is why a spec is free to stay
  loose. Append-only, never edited or re-dated: per-project `references.tsv`
  rows pin these **by sha256**, so regenerating one breaks a committed
  reproducibility pin.

## Rebuild commands

> **Note:** `bio.yml` no longer solves on linux01 — see the caveat under "pip
> installs on this box" below, and use the `--explicit` cache route instead.
> Specs (`sc.yml`, `meth.yml`) do still solve here.

Current env (matches what's live now):
```
conda env create -n bio -f bio.yml
```

A historical env (for reproducing an old project — look up which lock the
project pinned in its `references.tsv` `bio_env_lock` row):
```
conda env create -n bio_old -f lock/bio.2026-04-28.yml
```

## When to update

- **An export** (`bio`, `rna`, `primer`, `pwm2`, `pwm2-pb`) — re-export every
  time you `conda install` / `conda update` / `pip install` into it:
  ```
  conda env export -n bio > scripts/env/bio.yml
  ```
- **A spec** (`sc`, `meth`) — hand-edit it instead, adding the package under the
  right comment group. `check_env_drift.sh` flags the env once its `conda-meta`
  is newer than the yml, and editing the file is what clears the flag; committing
  a lock alone will not.
- **Snapshot a new lock** before any major upgrade (e.g. bumping HOMER,
  MACS3, bowtie2):
  ```
  cp scripts/env/bio.yml scripts/env/lock/bio.$(date +%F).yml
  ```
  Then commit both files together.

## pip installs on this box (glibc 2.23) — read before `pip install`

linux01 runs **glibc 2.23**. Most modern manylinux wheels are built against
glibc 2.25/2.28 and fail at import with `version 'GLIBC_2.xx' not found`.
Conda-forge builds link the env's own compatible libs and are always preferred.

**A failed pip-wheel import is NOT harmless.** It can leave the dynamic linker
holding a half-loaded dependency, and an unrelated call minutes later then
segfaults inside `ld-2.23.so`. Real case (diagnosed 2026-08-26): a pip
`pyarrow 24.0.0` wheel in `bio` could not import, and that alone caused
intermittent `segfault at b0 ... in ld-2.23.so` crashes in deepTools
(`plotHeatmap`) and plain pandas/numpy scripts from 2026-06 onward — plus it
silently broke `import sklearn` in `bio` for three months, because
`sklearn/utils/fixes.py` guards its `import pyarrow` with `except
ModuleNotFoundError` and the GLIBC failure raises a plain `ImportError` that
escapes the guard. Fix was `pip uninstall pyarrow`; nothing in `bio` needed it.

**Rules:**

1. Prefer `mamba install -c conda-forge <pkg>`. Only use pip when there is no
   conda build (pure-python packages are fine).
2. **Import-test immediately after any pip install**, in a fresh process:
   `$HOME/miniforge3/envs/<env>/bin/python -c "import <pkg>"`. A wheel that
   installs cleanly can still be unusable here.
3. Declare it afterwards (see "When to update" above): re-export an export,
   hand-edit a spec. `check_env_drift.sh` reports pip packages installed but
   absent from an **export**-style yml (`[pip-drift]`), and any pip package
   shipping compiled extensions that cannot import (`[pip-BROKEN]`).
4. Never set `CONDA_OVERRIDE_GLIBC=2.28` to force a package in. It silences the
   install-time solver check only; the ELF symbol requirement is still enforced
   by `ld.so` at runtime, which is exactly how these failures happen.

**Known recurrence vector:** MultiQC's *PyPI* metadata carries an unconditional
`Requires-Dist: pyarrow`, so `pip install -U multiqc` inside `bio` will silently
re-pull the broken wheel and re-break sklearn. The *bioconda* package correctly
depends on `polars-lts-cpu` instead — upgrade MultiQC with
`mamba install -n bio -c bioconda multiqc`, never with pip.

**Caveat on the rebuild commands above:** `conda env create -n bio -f bio.yml`
no longer solves on this box (tested 2026-07-24) — `r-base`/`pango` now resolve
to builds requiring glibc 2.28, and the glibc-2.23-compatible builds the live
env uses have aged out of channel repodata. Maintain `bio` **in place** with
targeted installs. To clone it exactly, bypass the solver via the local package
cache: `conda list -n bio --explicit > explicit.txt` then
`conda create -n <new> --file explicit.txt` (conda-only; re-add pip packages
separately).

Locks are append-only — never delete an old one. Old projects reference
them by filename.
