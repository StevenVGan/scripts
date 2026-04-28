# `scripts/env/` — conda environment tracking

The `bio` conda env is the single environment all pipelines (cutrun, csRNA,
proseq, atac) and tools rely on. This directory tracks its contents over
time so old projects remain reproducible after upgrades.

## Files

- **`bio.yml`** — current `bio` env, exported via `conda env export`.
  Re-export whenever you install/upgrade tools.
- **`lock/bio.YYYY-MM-DD.yml`** — dated snapshots. Save one before any
  major upgrade. Each per-project `references.tsv` records which lock
  file was active when the project's pipeline ran.

## Rebuild commands

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

- **Re-export `bio.yml`** every time you `conda install` / `conda update` /
  `pip install` into the `bio` env:
  ```
  conda env export -n bio > scripts/env/bio.yml
  ```
- **Snapshot a new lock** before any major upgrade (e.g. bumping HOMER,
  MACS3, bowtie2):
  ```
  cp scripts/env/bio.yml scripts/env/lock/bio.$(date +%F).yml
  ```
  Then commit both files together.

Locks are append-only — never delete an old one. Old projects reference
them by filename.
