# `~/work/` conventions

Single source of truth for layout, gitignore rules, file templates, and
workflows in `~/work/`. **AI agents and human collaborators: read this doc
before answering "where does X go?" or "how do I structure Y?"**

`~/work/CLAUDE.md` describes what `~/work/` *is* (a wet-lab analysis
filesystem tree, not a single repo). This doc describes how it's
*organized*.

---

## §1 — Quick reference (decision table)

| If you want to… | Do this |
|---|---|
| Set up a new sequencing project | See §setup. Copy `scripts/pipeline/<assay>/` into `seq/<assay>/<project>/script/pipeline/`, write `samples.tsv` (§4), `references.tsv` is auto-emitted by `5_qc.sh`, follow per-project README template (§4). |
| Add a one-off analysis script to project P | Put it in `seq/.../P/script/local/` and commit. |
| Same script needed in 2+ projects | Promote to `scripts/pipeline/tools/`. See §5. |
| New analysis spanning ≥2 projects | New repo under `seq/_joint/<name>/`. See §6. |
| Add a figure to a project repo | Drop in `figures/`, add a row to `figures/manifest.md`. Milestone-only — see §gitignore-rules. |
| Snapshot a final MultiQC report | Copy `multiqc/multiqc_report.html` → `figures/multiqc_YYYY-MM-DD.html`, add manifest row. |
| Track which env / refs were used | `scripts/env/lock/` for env snapshots; per-project `references.tsv` (auto-emitted) records which were active. |
| Use PWM motif analysis outputs in a seq project | Add a `dependencies.tsv` row in the seq project. See §7. |
| Upgrade the `bio` conda env | Re-export `bio.yml` and snapshot a new `lock/bio.<date>.yml`. See §8. |
| Run an existing project's pipeline | `cd seq/.../<project>/script && ./run_all.sh`. Toggles via env (e.g. `RUN_TRIM=0 ./run_all.sh`). |
| Make a change to a pipeline (`scripts/pipeline/<assay>/`) | Edit the source here. Per-project copies are stale by design until they're refreshed. |
| Find what's uncommitted across all repos | `~/work/gitall.sh` (when adopted) or `find ~/work -name .git -prune` + manual `git status`. |

---

## §2 — Repo layout map

```
~/work/
├── scripts/                # GIT REPO (source-of-truth pipelines + tools)
│   ├── env/                # conda env tracking (bio.yml, lock/bio.<date>.yml)
│   ├── CONVENTIONS.md      # this file
│   ├── pipeline/
│   │   ├── cutrun/         # standard CUT&RUN: 0_config.sh + numbered steps
│   │   ├── csRNA/          # csRNA fork (post-trim MultiQC gate, strand)
│   │   ├── proseq/         # PRO-seq fork (strand flip, pausing index)
│   │   ├── atac/           # ATAC-seq pipeline
│   │   └── tools/          # standalone utilities (heatmap.sh, peak_ops.sh,
│   │                       # go_enrichr.py, annotation_pie.py, etc.)
│   ├── experimental/       # WIP — peak-caller comparisons etc.
│   ├── legacy/             # older MACS2 single-file pipeline; reference only
│   └── project_archive/    # one-off scripts from finished projects
│
├── seq/<assay>/<project>/  # GIT REPO per active project
│   ├── README.md           # follows template (§4)
│   ├── 0_config.sh         # tracked: project-tuned
│   ├── peakcall_groups.tsv # tracked: ip/control/name/type pairing
│   ├── samples.tsv         # tracked: sample sheet (§4)
│   ├── references.tsv      # tracked: auto-emitted by 5_qc.sh (§4)
│   ├── script/
│   │   ├── pipeline/       # COPIES of scripts/pipeline/<assay>/ — NOT tracked
│   │   ├── local/          # tracked: project-specific scripts (incl link_fastq MAP_FILE)
│   │   └── run_all.sh      # tracked: project-tweaked driver
│   ├── analysis/<name>/    # tracked: within-project deeper analyses
│   ├── qc/summary.tsv      # tracked (when emitted)
│   ├── figures/
│   │   ├── manifest.md     # tracked
│   │   └── *.pdf|*.png|multiqc_<date>.html  # tracked: MILESTONE only
│   ├── data/, cleandata/, align/, peaks/, multiqc/, logs/, scratch/   # IGNORED
│
├── seq/_joint/<name>/      # GIT REPO per joint analysis (any source mix)
│   ├── README.md           # research question + which projects feed in
│   ├── sources.tsv         # tracked: source projects + pinned SHAs (§6)
│   ├── script/local/
│   ├── inputs/             # symlinks into source projects' outputs
│   ├── figures/manifest.md
│   ├── data/, scratch/     # IGNORED
│
├── PWM motif analysis/     # EXISTING independent git repo (notebooks/PWM tooling)
│
├── ref/                    # NOT tracked: genomes, blacklists, indexes (huge)
└── raw_seq/                # NOT tracked: raw FASTQs from IGM/ENA
```

---

## §3 — Gitignore rules

**Always ignored** (across all project / joint repos):

```
data/
cleandata/
align/
peaks/
multiqc/
logs/
scratch/
*.bam
*.bai
*.bigwig
*.bw
*.fastq.gz
*.fq.gz
script/pipeline/      # truth lives in scripts/ repo
.DS_Store
*.swp
__pycache__/
```

**Always tracked** (do not ignore):
- `0_config.sh`, `peakcall_groups.tsv`, `samples.tsv`, `references.tsv`,
  `README.md`, `script/local/*`, `script/run_all.sh`, `analysis/`,
  `qc/summary.tsv`, `figures/manifest.md`, `figures/*.pdf|*.png|*.html`.

**Rationale:**
- Code, config, and small text manifests → tracked. They're the
  reproducibility surface.
- Sequencing data and intermediates → ignored. They're huge and
  regenerable from raw FASTQ + tracked code + the env lock.
- `scratch/` → ignored. It's the exploratory dumping ground; if something
  there proves out, it graduates to `script/local/` and `figures/`.
- `script/pipeline/` → ignored at the project level; the source of truth
  is the `scripts/` repo. Per-project copies are pinned by `.pipeline_version`
  (when adopted) or recorded in `references.tsv`.

---

## §4 — File templates

### `samples.tsv` (per project, tracked)

The sample sheet — one row per sample. Subsumes any `link_fastq`
`MAP_FILE` so data provenance lives with sample metadata.

```
sample_id            target    condition    replicate    raw_fastq_R1                          raw_fastq_R2                          sequencing_run    notes
MCF7_ERa_DMSO_rep1   ERa       DMSO         1            raw_seq/260401_IGM/SG13_S1_R1.fq.gz   raw_seq/260401_IGM/SG13_S1_R2.fq.gz   260401_IGM        Priyanka submission
MCF7_ERa_E2_rep1     ERa       E2_1h        1            raw_seq/260401_IGM/SG14_S2_R1.fq.gz   raw_seq/260401_IGM/SG14_S2_R2.fq.gz   260401_IGM        Priyanka submission
```

Use `-` for `raw_fastq_R2` on single-end and for any column that doesn't
apply. Use `notes` to flag pre-publication or external-collaborator data
(prefix `EMBARGO:` so a pre-push scan can catch it).

### `references.tsv` (per project, tracked, auto-emitted)

Records which reference files (and which env lock) the pipeline run used.
Auto-emitted at the end of `5_qc.sh`. You don't author this by hand.

```
key             path                                                                       size_bytes    mtime_iso             sha256_or_dash
genome_fasta    /mnt/share/archive/bkup/ref/align/bowtie2/hg38_noalt/hg38.fa               3209286105    2024-08-12T10:13:42   -
bowtie2_index   /mnt/share/archive/bkup/ref/align/bowtie2/hg38_noalt/hg38                  -             -                     -
blacklist       /mnt/home/digan/work/ref/blacklist/hg38/hg38-blacklist.v2.bed              1234567       2023-04-01T00:00:00   abc123def456…
bio_env_lock    /mnt/home/digan/work/scripts/env/lock/bio.2026-04-28.yml                    23048         2026-04-28T05:08:00   -
```

`sha256` only for small files (blacklists, custom BEDs). Big files get
size+mtime — catches "someone replaced the file" without minutes of hashing.

### `sources.tsv` (per joint repo, tracked)

Records source projects + pinned SHAs. Authored manually when the joint
is created; refresh SHAs when you re-pull from a source.

```
name              path                                              git_sha    role
ERa_OGG1_KD       ../../CUTRUN/260115_CnR_ERa_OGG1_MCF7_KD_Priyanka  abc1234   knockdown
ERa_OGG1_nonKD    ../../CUTRUN/260401_CnR_ERa_OGG1_MCF7_Priyanka     def5678   control
```

### `dependencies.tsv` (per project consuming PWM motif analysis or other external repo)

```
name                  path                              git_sha    used_for
pwm_motif_analysis    ../../PWM motif analysis           a1b2c3d   cleavage site motif scoring
```

### `figures/manifest.md`

```
# Milestone figures

| File | Generated by | Date | Notes |
|---|---|---|---|
| er_p65_cobinding_heatmap.pdf | script/local/cobinding_heatmap.sh | 2026-04-15 | for Apr lab meeting |
| multiqc_2026-04-28.html      | scripts/pipeline/cutrun/5_qc.sh    | 2026-04-28 | post-restructure final QC |
```

### Per-project `README.md`

```markdown
# <project-name>

**Assay:** <CUTRUN | ChIPseq | PROseq | csRNA | CnT | ATAC | …>
**Genome:** <hg38 | hg19 | mm10 | …>
**SE/PE:** <SE | PE>
**Owner / submitter:** <name>
**Date received / batch:** <YYYY-MM-DD or batch id>
**Status:** <active | analysis-complete | published | paused>
**Upstream source:** <GEO accession / IGM submission / collaborator>

## What this project is
<2–3 sentences: biological question, IP target(s), comparison being made>

## Key analyses
- `analysis/<name>/` — <one-line description>
- joint repos this project feeds: <`seq/_joint/<name>/` or "none">
- joint repos this project consumes: <list or "none">

## Run
- Pipeline: `scripts/pipeline/<assay>/` (see `references.tsv` `bio_env_lock`)
- Standard invocation: `cd script && ./run_all.sh`
```

---

## §5 — Promotion path (project local → shared tools)

When a `script/local/` script gets copy-pasted into a 2nd project (or is
clearly general-purpose), promote it.

```
project script/local/  →  scripts/pipeline/tools/   (or scripts/pipeline/<assay>/ if it's a numbered pipeline step)
```

**Mechanics:**
1. In `scripts/`: `git mv` (or copy + delete) the script to its new home.
   Adjust shebang/paths if needed. Commit with message
   `promote <name> from <originating_project>`.
2. In the originating project's `script/local/`: replace with a one-line
   wrapper that `exec`s the promoted version, **or** delete and update
   callers to use the shared path. The wrapper pattern (used by existing
   `link_fastq.sh` per-project wrappers) keeps invocation paths stable.

**Topic subfolders inside `tools/`** are encouraged when a methodology
involves multiple related scripts that belong together. Examples:

```
scripts/pipeline/tools/
├── cleavage/         # cleavage-site methodology (BAM→cut bigwig, motif cutoff sweeps, footprint plots)
├── cobinding/        # ChIP/CUT&RUN cobinding consensus + downstream
├── heatmap.sh        # generic deepTools heatmap (flat tool — fine to leave at top)
├── go_enrichr.py     # generic Enrichr GO (flat tool)
└── ...
```

Create a topic subfolder the moment you'd otherwise be putting ≥2 related
scripts side-by-side at the top of `tools/`. Keep flat single-script
utilities at the top level.

The `scripts/experimental/` middle tier exists for things that are
"promoted but not yet stable enough for `tools/`" — use sparingly. Most
promotions go straight to `tools/` for solo work.

---

## §6 — Joint analyses

When an analysis depends heavily on data from ≥2 projects, it gets its
own home: `seq/_joint/<name>/`. The leading `_` sorts joint dirs apart
from per-project dirs in `ls`.

**Naming — recurring analysis families share a prefix.** Lead with the
analysis subject (e.g. `cleavage_sites_*`, `cobinding_*`,
`motif_density_*`), then the dataset combo or biological focus
(`_ERa_OGG1_KD_combined`, `_MCF7_ER_p65`, etc.). Why:

- Related joints cluster alphabetically — `ls seq/_joint/cleavage_sites_*`
  shows the whole family.
- Greppable across `sources.tsv` files when asking "what other joints
  consume project X?"
- Triggers the natural promotion path: when ≥2 joints in a family share
  scripts, those scripts belong in `scripts/pipeline/tools/<family>/`
  (see §5 topic subfolders).

**One convention covers all joints**, regardless of whether sources are
same-assay or cross-assay. (No separate `_integrated/` directory.) The
README documents what assays are involved.

**Setup:**
1. `mkdir -p seq/_joint/<name>/{script/local,inputs,figures,data,scratch}`
2. `git init` in the new dir.
3. Drop in the standard project `.gitignore` (§3).
4. Author `sources.tsv` (§4) listing source projects + their current SHAs.
5. Populate `inputs/` with symlinks into source projects' outputs:
   `ln -s ../../<assay>/<source-project>/peaks/macs3/<file>.bed inputs/<role>_peaks.bed`
6. Author `README.md` (research question, what's inside).
7. Move analysis scripts into `script/local/`.
8. First commit, push to GitHub.

**When sources update:** bump SHAs in `sources.tsv`, decide whether to
re-run analysis, commit. (No automation today; manual is fine for now.)

---

## §7 — PWM motif analysis usage

`PWM motif analysis/` is an independent git repo
(github.com/StevenVGan/PWM-motif-analysis). It is *not* a sequencing
project and does not consume the standard pipelines. Several seq projects
consume its outputs or scripts.

**To use PWM in a seq project:** add `dependencies.tsv` (§4) recording the
PWM commit SHA you built against. Update the SHA whenever you re-pull from
a newer PWM commit.

**If a PWM utility matures into a generally-useful seq tool:** promote it
into `scripts/pipeline/tools/` per §5. Commit message:
`promote <name> from PWM motif analysis @ <sha>`.

---

## §8 — Conda env (`bio`)

The `bio` env has every tool the pipelines invoke (bowtie2, MACS3, HOMER,
deepTools, samtools, multiqc, R packages, …). Tracked in
`scripts/env/`. Full instructions in `scripts/env/README.md`.

**Quick reference:**
- Current: `scripts/env/bio.yml`
- Snapshots: `scripts/env/lock/bio.YYYY-MM-DD.yml`
- Re-export after install/upgrade:
  `conda env export -n bio > scripts/env/bio.yml`
- Snapshot before major upgrade:
  `cp scripts/env/bio.yml scripts/env/lock/bio.$(date +%F).yml`
- Per-project pin: each project's `references.tsv` `bio_env_lock` row
  records which lock was active during its pipeline run.

Locks are append-only. Never delete old ones.

---

## §9 — Reference data (`ref/`)

`~/work/ref/` holds genomes, bowtie2 indexes, blacklists, ER enhancer
BEDs, etc. Too big to track in git. Symlinked / shared via filesystem.

Per-project `references.tsv` (§4, auto-emitted by `5_qc.sh`) records the
*paths*, *sizes*, *mtimes*, and (for small files) *sha256* of the
references that were active when the pipeline ran. This is the drift
detector: if a blacklist file gets replaced silently, the next
`references.tsv` regeneration will show a different sha256 / mtime.

**Bowtie2 index path convention:**
`/mnt/share/archive/bkup/ref/align/bowtie2/<genome>_noalt/<genome>`

---

## §10 — Logs

Pipeline logs go to `<project>/logs/<step>_<timestamp>.log` via `log_start`
in `0_config.sh`. The `logs/` dir is always gitignored.

**Auto-prune:** `log_start` keeps the **3 most recent** logs per step
(`LOG_KEEP_N=3`, override per-invocation if needed). Older logs auto-delete.
Tight cap because pipeline reruns are common during debugging.

**Exempt from auto-prune:** anything under `<project>/logs/keep/`. Move a
log there manually if you want to preserve it (failed-run diagnosis,
benchmark, etc.). Auto-prune does not touch `logs/keep/`.

**Exploration vs. canonical runs:** `<project>/scratch/` is gitignored and
*not* subject to log auto-prune. Send exploratory commands' stdout/stderr
there so the 3-most-recent cap on `logs/` only counts canonical pipeline
runs.

---

## §setup — Setting up a new sequencing project

1. `mkdir -p seq/<assay>/<new_project>/{data,script/local,figures,qc,analysis}`
2. Copy pipeline: `cp scripts/pipeline/<assay>/* seq/<assay>/<new_project>/script/pipeline/`
   (Or: once `setup_project.sh` exists, just run it.)
3. Edit `seq/<assay>/<new_project>/script/0_config.sh`:
   - Set `BASE` to the project root.
   - Set `GENOME`, `SE` if non-default.
4. Author `samples.tsv` (template §4).
5. Author `peakcall_groups.tsv` (ip/control/name/type rows).
6. Author `README.md` (template §4).
7. Author `link_fastq.sh` MAP_FILE under `script/local/` (or fold mapping
   into `samples.tsv`). Run `script/local/link_fastq.sh` to populate `data/`.
8. `git init`, drop in the standard `.gitignore` (§3), first commit.
9. Run pipeline: `cd script && ./run_all.sh`. `references.tsv` is
   auto-emitted by `5_qc.sh`.
10. Add remote, push.
