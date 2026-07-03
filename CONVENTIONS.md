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
| Set up a new sequencing project | Run `scripts/setup/new_project.sh <pipeline> <name>` (§setup) — it stamps the flat `script/`, `.gitignore`, and template TSVs. Then fill `samples.tsv` / `peakcall_groups.tsv` (§4); `references.tsv` auto-emits from `5_qc.sh`. |
| Add a one-off analysis script to project P | Put it in `seq/.../P/script/local/` and commit. |
| Same script needed in 2+ projects | Promote to `scripts/pipeline/tools/`. See §5. |
| New analysis spanning ≥2 projects | New repo under `seq/_joint/<name>/`. See §6. |
| Add a wet-lab experiment derived from a project (qPCR panel, knockdown, validation) | New dir under `seq/.../P/experiments/<name>/`. Wet-lab work generates *new* data — keep it out of `analysis/` (which is in-silico re-analysis of the deposit). See §experiments. |
| Start a multi-step in-silico analysis of a project | New dir under `seq/.../P/analysis/<topic>_v1/` with numbered Python steps. See §11. |
| Materially change an analysis approach (feature set, model, label space) | Bump version: copy `<topic>_v1/` → `<topic>_v2/`, add `Supersedes:` line to v2's README and `Superseded by:` to v1's. See §11. |
| Hand an artifact from one analysis to another | Pin the producer path in the consumer's `00_config.sh` via env var (e.g. `FOOTPRINT_TSV=...`). See §11. |
| Set up a multiome or scRNA project | Different toolchain (scanpy/SnapATAC2 in the `sc` env, not `bio`). See §12. |
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
│   ├── setup/              # new_project.sh scaffolder + PROVISION.md (2nd-node bootstrap)
│   ├── CONVENTIONS.md      # this file
│   ├── pipeline/
│   │   ├── cnr/            # standard CUT&RUN: 0_config.sh + numbered steps
│   │   ├── csrna/          # csRNA fork (post-trim MultiQC gate, strand)
│   │   ├── pro/            # PRO-seq fork (strand flip, pausing index)
│   │   ├── atac/           # ATAC-seq pipeline
│   │   ├── multiome/        # single-cell / multiome scaffold
│   │   ├── templates/      # canonical scaffolding templates (gitignore, TSVs, README)
│   │   └── tools/          # standalone utilities (heatmap.sh, peak_ops.sh,
│   │                       # go_enrichr.py, annotation_pie.py, etc.)
│   ├── experimental/       # WIP — peak-caller comparisons etc.
│   ├── legacy/             # older MACS2 single-file pipeline; reference only
│   └── project_archive/    # one-off scripts from finished projects
│
├── seq/<assay>/<project>/  # GIT REPO per active project
│   ├── README.md           # follows template (§4)
│   ├── peakcall_groups.tsv # tracked: ip/control/name/type pairing
│   ├── samples.tsv         # tracked: sample sheet (§4)
│   ├── references.tsv      # tracked: auto-emitted by 5_qc.sh (§4)
│   ├── link_sample.tsv     # tracked: link_fastq MAP_FILE — data provenance
│   ├── multiqc_config.yaml # tracked: project-specific MultiQC config
│   ├── script/             # FLAT (no script/pipeline/ subfolder — see §3 lean note)
│   │   ├── 0_config.sh     # tracked: project-tuned
│   │   ├── run_all.sh      # tracked: project-tweaked driver
│   │   ├── link_fastq.sh   # tracked: project wrapper around scripts/pipeline/tools/prep/link_fastq.sh
│   │   ├── 1_trim_qc.sh    # NOT tracked — copy of scripts/pipeline/<assay>/1_trim_qc.sh
│   │   ├── 2_bowtie2.sh    # NOT tracked — copy
│   │   ├── 3_homer_tags.sh # NOT tracked — copy
│   │   ├── 4.1_*.sh, 4.2_*.sh, 4.3_*.sh   # NOT tracked — copies
│   │   ├── 5_qc.sh         # NOT tracked — copy
│   │   └── local/          # tracked: project-specific bespoke scripts
│   ├── analysis/<name>/    # tracked: within-project deeper analyses of the deposit (heavy outputs ignored — see §3)
│   ├── experiments/<name>/ # tracked: wet-lab experiments derived from the analysis (qPCR/knockdown/validation) — panel/protocol/small results; raw data ignored (see §3, §experiments)
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

> **Lean-track note (Phase C, 2026-04-28):** The original plan called for
> moving numbered step scripts into `script/pipeline/` (a subfolder) and
> gitignoring the subfolder. In practice we kept `script/` flat and use
> a glob rule to ignore the numbered scripts in place. Same principle
> ("pipeline copies aren't tracked"), zero edits to sourcing paths,
> smaller per-project diff. All 12 active projects use this layout.

**Standard `.gitignore` (per project repo)** — canonical source is
`scripts/pipeline/templates/gitignore` (stamped by `new_project.sh`).
**One pattern per line**: git treats a line with unescaped spaces as a *single*
literal pattern, so multi-pattern lines silently match nothing.

```gitignore
# Data + intermediates (regenerable from raw FASTQ + tracked code + bio env lock)
data/
cleandata/
align/
peaks/
multiqc/
logs/
scratch/

# Common project-specific intermediate dirs (varies; extend per project)
plot/
old_heatmap/
oridata/
pausing/
tss/
analysis/inputs/

# Heavy intermediate dirs INSIDE analysis/ at any depth (regenerable)
analysis/**/bigwig/
analysis/**/matrices/
analysis/**/motif/
analysis/**/qc/
analysis/**/seq/
analysis/**/tables/
analysis/**/regions/
analysis/**/bed/
analysis/**/tracks_mean/
analysis/**/figures/
analysis/**/macs/
analysis/**/homer/
analysis/**/data/
analysis/**/heatmap/
analysis/**/merged_bw/
analysis/**/bw_farm/
analysis/**/peaks/
analysis/**/logs/
analysis/**/classify/

# Heavy raw data inside experiments/ (qPCR exports, instrument files, gel/microscopy images)
experiments/**/data/
experiments/**/raw/
experiments/**/logs/

# HOMER tag directories (large per-chr .tags.tsv files; regenerable)
**/*_tagdir/
**/pooled_*_tagdir/
**/tagdirs/
*.tags.tsv

# Binary outputs (anywhere) — ONE PER LINE
*.bam
*.bai
*.bigwig
*.bw
*.bedgraph
*.bg
*.fastq.gz
*.fq.gz
*.npz

# Compressed matrices
*.mat.gz
*.matrix.gz
matrix_*.gz

# HOMER motif outputs (numerous, regenerable)
*.motif
*.motifs
*.motifs8
*.motifs10
*.motifs12
*.svg
*.html
**/homerResults/
**/knownResults/

# MACS pileup intermediates
*_treat_pileup.bdg
*_control_lambda.bdg
*_peaks.xls

# Loose figure outputs scattered at subdir top-level
analysis/*/figure_*.png
analysis/*/figure_*.pdf
analysis/*/*_heatmap.png
analysis/*/*_profile_*.png

# Caches + per-analysis run logs
**/bed_cache/
analysis/**/*.log
*.log.tsv

# Pipeline copies (truth in scripts/pipeline/<assay>/ — scripts/ repo)
script/[1-9]*.sh
script/[1-9]*.py
script/README.md
script/.lintr

# Editor / OS / Python — ONE PER LINE
.DS_Store
*.swp
__pycache__/
*.pyc
.ipynb_checkpoints/
```

**Always tracked** (do not ignore):
- `0_config.sh`, `run_all.sh`, `link_fastq.sh` (project wrappers)
- `peakcall_groups.tsv`, `samples.tsv`, `references.tsv`, `link_sample.tsv`,
  `multiqc_config.yaml` (project metadata)
- `README.md`
- `script/local/*` (project-specific bespoke scripts)
- `analysis/` script files (.sh, .py, .R) and small-text outputs
- `qc/summary.tsv` (when emitted)
- `figures/manifest.md` and milestone figures (PDF/PNG/MultiQC HTML snapshots)

**Rationale:**
- Code, config, and small text manifests → tracked. They're the
  reproducibility surface.
- Sequencing data and intermediates → ignored. They're huge and
  regenerable from raw FASTQ + tracked code + the env lock.
- `scratch/` → ignored. It's the exploratory dumping ground; if something
  there proves out, it graduates to `script/local/` and `figures/` with a
  manifest entry.
- `script/[1-9]*.sh` → ignored at the project level; the source of truth
  is the `scripts/` repo. Per-project copies match upstream by convention.
- `analysis/**/<heavy-dir>/` → ignored to prevent gigabytes of derived
  outputs (HOMER tagdirs, matrices, motif/, bigwigs, etc.) from bloating
  per-project repos.

**Per-project extras** are common — e.g. `analysis/**/dinucleotide_conservation_long.tsv`,
`analysis/**/refFlat*.txt.gz`, etc. Add these to a project's `.gitignore`
when the standard rules let through unwanted files. Audit with
`git ls-files | xargs ls -la | sort -k5n -r | head` after first commit.

**Non-standard layouts (rare):** project #9 (`240216_240701_csRNA`) uses
`scripts/csRNA/` instead of `script/`. Its gitignore reflects this with
`scripts/csRNA/[1-9]*.sh`. Document any such deviation in the project's
`README.md`.

---

## §4 — File templates

### `samples.tsv` (per project, tracked)

The sample sheet — one row per sample; the human-authored record of sample
metadata + raw-FASTQ provenance. Distinct from `link_sample.tsv` (below), the
terse `prefix→newname` map that `link_fastq` actually consumes.

```
sample_id            target    condition    replicate    raw_fastq_R1                          raw_fastq_R2                          sequencing_run    notes
MCF7_ERa_DMSO_rep1   ERa       DMSO         1            raw_seq/260401_IGM/SG13_S1_R1.fq.gz   raw_seq/260401_IGM/SG13_S1_R2.fq.gz   260401_IGM        Priyanka submission
MCF7_ERa_E2_rep1     ERa       E2_1h        1            raw_seq/260401_IGM/SG14_S2_R1.fq.gz   raw_seq/260401_IGM/SG14_S2_R2.fq.gz   260401_IGM        Priyanka submission
```

Use `-` for `raw_fastq_R2` on single-end and for any column that doesn't
apply. Use `notes` to flag pre-publication or external-collaborator data
(prefix `EMBARGO:` so a pre-push scan can catch it).

### `peakcall_groups.tsv` (per project, tracked)

Peak-calling jobs, one per row. Read by `4.1_peak_macs3.sh` + `4.2_peak_homer.sh`.
TAB-separated, **no header row**; blank and `#`-comment lines are skipped.

```
# ip_bam                      control_bam              name              type
SG01_ERa_E2_rep1_sorted.bam   SG02_IgG_E2_sorted.bam   SG01_ERa_E2_rep1  TF
SG05_H3K27ac_rep1_sorted.bam  -                        H3K27ac_rep1      Histone
```

`control_bam` = `-`/`none`/`NA` for no control; `type` = `TF` (narrow) | `Histone`
(broad). The csRNA fork additionally uses `tss_groups.tsv` (same shape minus the
`type` column), read by `4.3_tss_csrna.sh` → `findcsRNATSS.pl`; the PRO fork
replaces that step with `4.3_pausing_divergent.sh` and uses no groups file.

### `link_sample.tsv` (per **project root**, tracked)

The `link_fastq` `MAP_FILE`: raw FASTQ prefix → canonical sample name. Canonical
location is the project root. TAB-separated, no header.

```
SG_273   SG273_CnR_ERa_MCF7_Ctrl_E2_rep1
SG_274   SG274_CnR_ERa_MCF7_Ctrl_E2_rep2
```

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

`sha256` only for files **< 10 MB** (blacklists, custom BEDs). Larger files get
size+mtime — catches "someone replaced the file" without minutes of hashing.

### `sources.tsv` (per joint repo, tracked)

Records source projects + pinned SHAs. Authored manually when the joint
is created; refresh SHAs when you re-pull from a source.

```
name              path                                              git_sha   role
ERa_OGG1_KD       ../../cnr/CnR_260115_ERa_MCF7_OGG1-KD_Priyanka    abc1234   knockdown
ERa_OGG1_nonKD    ../../cnr/CnR_260401_ERa_MCF7_OGG1-inhi_Priyanka  def5678   control
```

### `dependencies.tsv` (cross-project / external-artifact ledger, tracked)

The workspace dependency graph: any time a project or analysis consumes another
project's artifact (an external repo like PWM, a produced h5ad / DEG table /
signature, or a shared helper), record it — one row per artifact, fixed schema
`name / path / git_sha / used_for`. **Fill `git_sha`** with the producer's commit
(use `-` only when the source isn't git-init'd, and say so in `used_for`).
sc/multiome analyses never run `5_qc.sh`, so they get no auto `references.tsv`;
give each one a `sc_env_lock` row pinning the `scripts/env/lock/sc.<date>.yml` it
ran under, so the analysis still records its environment.

Paths are relative to the project root where `dependencies.tsv` lives:
`~/work`-level repos (`PWM motif analysis`, `scripts`) are three `../` up; a
sibling seq project is `../../<assay>/...`.

```
name                  path                                                       git_sha   used_for
pwm_motif_analysis    ../../../PWM motif analysis                                6c75a25   cleavage-site motif scoring
shared_signature      ../../scrna/<project>/analysis/<name>_v1/sig.tsv           b7e912d   score cells against a shared signature
sc_env_lock           ../../../scripts/env/lock/sc.2026-05-30.yml               -         sc env lock (dated filename is the pin)
```

### `figures/manifest.md`

```
# Milestone figures

| File | Generated by | Date | Notes |
|---|---|---|---|
| er_p65_cobinding_heatmap.pdf | script/local/cobinding_heatmap.sh | 2026-04-15 | for Apr lab meeting |
| multiqc_2026-04-28.html      | scripts/pipeline/cnr/5_qc.sh    | 2026-04-28 | post-restructure final QC |
```

### Per-project `README.md`

```markdown
# <project-name>

**Assay:** <CUT&RUN | ChIPseq | PRO-seq | csRNA | CnT | ATAC | multiome | scRNA | …>
**Genome:** <hg38 | hg19 | mm10 | …>
**SE/PE:** <SE | PE>
**Owner / submitter:** <name>
**Date received / batch:** <YYYY-MM-DD or batch id>
**Status:** <active | analysis-complete | published | paused>
**Upstream source:** <GEO accession / IGM submission / collaborator>

## What this project is
<2–3 sentences: biological question, IP target(s), comparison being made>

## Key analyses
- `analysis/<name>_v<N>/` — <one-line description> [+ `Supersedes:`/`Superseded by:` if relevant]
- `experiments/<name>/` — <wet-lab experiment, if any>
- joint repos this project feeds: <`seq/_joint/<name>/` or "none">
- joint repos this project consumes: <list or "none">

## Run
- Pipeline: `scripts/pipeline/<assay>/` (see `references.tsv` `bio_env_lock`) — bulk-seq only
- Standard invocation: `cd script && ./run_all.sh`
- Multiome / scRNA: no project-level pipeline. Analyses are self-driving under `analysis/<name>/script/run_all.sh`. See §11.
```

### Per-analysis `README.md` (inside `analysis/<name>/`)

```markdown
# <topic>_v<N>

**Supersedes:** [v<N-1>](../<topic>_v<N-1>/)         <!-- omit if v1 -->
**Superseded by:** [v<N+1>](../<topic>_v<N+1>/)       <!-- add when bumped -->
**Status:** <pending | running | complete | abandoned>
**Headline result:** <one sentence; numbers if you have them>

## Pipeline summary
<one-paragraph "binarize → TF-IDF → SVD → GMM" style description; the
authoritative source of truth is the numbered scripts in `script/`>

## Inputs
- `<env-var>=<default-path>` — <what it points at, which producer version>

## Outputs
- `processed/<key-table>.tsv` — <consumed by which downstream analysis>
- `figures/<key-figure>.pdf` — <appears in REPORT.md / manifest>
```

### `RESUME.md` (per-analysis, when handing off mid-state)

```markdown
# Resume state — <topic>_v<N>

**Last touched:** <YYYY-MM-DD>
**Currently in:** <step NN — what's running / blocked / waiting on>
**Next action:** <single concrete step the next session should take>

## What's done
- [x] 01–<NN-1> — <one-line summary>

## What's running
- step <NN> launched at <time> as `nohup … &`, PID <pid>, log `logs/<step>_<ts>.log`

## What's blocked / decisions pending
- <e.g. "the v3.3 cohort split dropped 14% of cells at QC — need Steven's call on filter strictness">
```

### `CHANGELOG.md` (per-analysis, append-only)

```markdown
# Changelog — <topic>

## v<N> (YYYY-MM-DD)
- <one-line summary of the change>
- **Supersedes** v<N-1> because <reason>

## v<N-1> (YYYY-MM-DD)
- …
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

### Repo-level graduation (rare, for tool families that outgrow `tools/`)

When a tool family grows into its own ontology + API + frozen release
bundle, it can graduate to its own **top-level repo** under `~/work/`,
alongside `PWM motif analysis/`. Trigger: a tool that started in
`analysis/<name>/` of a single project but is now generic across tissues
and consumed by ≥2 downstream analyses.

`tf_atlas/` is the exemplar: extracted May 2026 from a per-project
`analysis/tf_atlas/` when it ceased to be tissue-specific. Consumers import via

```python
sys.path.insert(0, '/mnt/home/digan/work/tf_atlas/script')
from atlas_api import TFAtlas
```

and depend on a frozen `RELEASE_v<N>.<M>/` bundle inside the tool repo so
consumers can pin a known-good snapshot independent of `main`. Layout:

```
~/work/<toolname>/
├── README.md
├── docs/RELEASE_v<N>_plan.md
├── script/
│   ├── atlas_api.py        # public entry point
│   ├── _utils.py           # shared helpers (see §11)
│   ├── pipeline/           # batch build steps
│   ├── orchestrate/        # wiring across phases
│   ├── api/                # F1/F2/F3-style query interfaces
│   ├── diagnostics/        # null-perm tests, regression suites
│   └── legacy/             # frozen old versions
├── data/atlas_artifacts/RELEASE_v<N>.<M>/   # frozen snapshot
└── tests / regression suite
```

A `RELEASE_v<N>.<M>/` bundle is **frozen and self-describing**: alongside the
data it carries its own `MANIFEST.tsv` (files + shas), `METHODS.md`,
`CHANGELOG.md`, and `NEGATIVE_RESULTS.md`, so it can be audited without the build
tree (tf_atlas's is the exemplar). **Consumers must pin to the bundle dir** and
record the tool's commit SHA in their `dependencies.tsv` — not float on live
`main` with only a soft version string, or the "known-good snapshot independent
of `main`" guarantee is lost.

Don't graduate prematurely. The bar is: ≥2 downstream consumers, an
internal versioning convention (§11) already needed, and a release-bundle
discipline mature enough to commit to.

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

## §experiments — Wet-lab experiments derived from a project

`analysis/` holds **in-silico re-analysis of a project's deposited data**.
When a finding motivates a **wet-lab experiment** (qPCR validation, knockdown,
luciferase, etc.), it generates *new* experimental data and belongs in a
sibling **`experiments/<name>/`** dir, not in `analysis/`. Keeping the two
apart preserves the meaning of `analysis/` (dry, regenerable from the deposit)
vs. `experiments/` (wet, new primary data).

**Layout:**

```
seq/<assay>/<project>/experiments/<name>/
├── README.md      # hypothesis, design, data provenance (which analysis/ dirs it derives from), caveats, status
├── panel.tsv      # tracked: gene/primer panel, protocol params, or design table
├── primers.tsv    # tracked: primer sequences (when designed)
├── results/       # tracked: small results tables, summary stats
├── figures/       # ignored by default (§3) — promote milestones like any project figure
├── data/, raw/    # IGNORED (§3) — qPCR exports, instrument files, images
```

**Conventions:**
- It lives **inside** the project whose analysis motivated it — provenance
  stays local (the README points at the `analysis/<name>/` tables it derives
  from via `../../analysis/...`).
- Same name-with-version discipline as analyses (`<topic>_v1`, `_v2`).
- Tracked: README, design/panel/primer/results text tables. Ignored: raw
  instrument data and heavy outputs (`experiments/**/{data,raw,logs}/` — §3).
- **Graduation:** if the wet-lab line grows its own sequencing (RNA-seq,
  CUT&RUN, etc.), it graduates to its own `seq/<assay>/<project>/`. A single
  qPCR/validation experiment does not — it stays a project-local experiment.

**Two homes, by provenance.** The above is for an experiment that derives from a
project's deposit — it lives *inside* that project. An experiment with **no**
sequencing deposit behind it (a standalone bench assay) instead goes in the
top-level **`~/work/experiments/`** repo as a dated sibling dir
(`<YYMMDD>_<slug>_v<N>/`), with shared qPCR helpers factored into
`experiments/_lib/` on their 2nd use. Routing test: "does this derive from a
deposit?" — yes → project-local; no → the top-level repo.

First instance (project-local): a `seq/<assay>/<project>/experiments/<name>_v1/`
qPCR validation panel derived from that project's `analysis/` finding.

---

## §7 — PWM motif analysis usage

`PWM motif analysis/` is an independent git repo
(github.com/StevenVGan/PWM-motif-analysis). It is *not* a sequencing
project and does not consume the standard pipelines. Several seq projects
consume its outputs or scripts.

**To use PWM in a seq project:** add a `dependencies.tsv` row (§4) recording the
PWM commit SHA you built against; update it when you re-pull. `dependencies.tsv`
is in fact the general ledger for *any* cross-project or external-repo
consumption (§4), not just PWM — one row per consumed artifact, `git_sha` filled.

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

**Preferred — one command:**

```bash
scripts/setup/new_project.sh <pipeline> <project_name>   # pipeline = cnr | csrna | pro | atac
#   --dest <subdir>   place under seq/<subdir>/  (default = pipeline; e.g. --dest cnt or --dest chip)
#   --genome hg19 --se --owner <name> --upstream <accession>
```

It stamps the full skeleton — the dir tree, the pipeline scripts (flat in
`script/`, git-ignored copies), the one-per-line `.gitignore`, the header
template TSVs, a token-filled `README.md`, and a `link_fastq.sh` that `exec`s
the shared prep script — then `git init` + stages (it does **not** commit).
`BASE` auto-derives from the project root, so `0_config.sh` needs no edit;
non-default `GENOME`/`SE` are baked in by the flags.

Then supply the biology the scaffolder can't:
1. `link_sample.tsv` — raw prefix → sample name; set `RAW_DIR` in `script/link_fastq.sh`.
2. `samples.tsv` — one row per sample (§4).
3. `peakcall_groups.tsv` — IP/control pairings + `TF`/`Histone` (§4).
4. `README.md` — fill "What this project is" + any remaining header TODOs.
5. `cd script && ./link_fastq.sh && ./run_all.sh` — `references.tsv` auto-emits from `5_qc.sh`.
6. Review (`git status`), commit, add remote (`gh repo create`), push.

---

## §11 — Per-project deep analyses (`analysis/<name>/`)

Per-project `analysis/<name>/` dirs hold *in-silico re-analyses of the
deposit* (vs `experiments/<name>/` which is wet-lab; §experiments). When
an analysis grows beyond a handful of scripts, it follows the same
config-driven + numbered-step discipline as the bash pipeline — but in
Python, not bash.

### Versioning

When an analysis approach changes materially (different feature set,
different statistical model, different label space), **bump the version
and keep the old one**:

```
analysis/
├── clustering_v3/            # canonical
├── clustering_v4/            # generalization test (kept; weaker but informative)
├── clustering_v5/            # ML feature ranking (kept; complementary)
├── clustering_v3.3_cohortB/  # cohort split of v3.3
└── de_analysis_v3/           # supersedes v2 (see README)
```

- `<topic>_vN` for major versions — a **new feature set, statistical model,
  label space, or input cohort** (the result could legitimately differ).
  `<topic>_vN_M` (`v3_5`, `v1.1_pooled`) for incremental tweaks that leave the
  method identity intact (a threshold, a seed, a plotting change).
- Suffix sub-variants with a context tag: `clustering_v3.3_cohortB/`
  vs `clustering_v3.3/` (cohort split), not `_v3.3a/_v3.3b/`.
- **Never delete a superseded version** without first folding a `REPORT.md`
  summary into the parent project's `CHANGELOG.md` or `docs/`. Old
  versions document the paths that didn't work.
- When v2 supersedes v1, add `**Supersedes:** [v1](../<topic>_v1/)` at
  the top of v2's README, and `**Superseded by:** [v2](../<topic>_v2/)`
  at the top of v1's README — visible without clicking through.
- **Retire, don't delete.** Move a superseded version into `analysis/_archive/`
  (the leading `_` sorts it aside, like `seq/_joint/`) and add a row to
  `analysis/_archive/ARCHIVE.md` naming the dir, what replaced it, and one line
  on why. Keeps the live `analysis/` listing readable while preserving the
  dead-end for the record; the `ARCHIVE.md` row also substitutes for a
  `Superseded by:` pointer when the retired dir has no README.

### Doc set inside `analysis/<name>/`

| File | When to use |
|---|---|
| `README.md` | Always. Headline result, status, pipeline summary, supersedes pointer. |
| `PLAN.md` | When the analysis spans ≥5 steps or multiple sessions. Living design doc. |
| `RESUME.md` | When handing off mid-state (own future-self after a break, or a subagent). |
| `REPORT.md` | Long-form writeup of the final result. May graduate to parent project's `docs/`. |
| `CHANGELOG.md` | Append-only history of version bumps and what changed. |
| `BACKLOG.md` | Open follow-ups that aren't urgent enough to block. |

Only `README.md` is required. Add the others when they earn their keep.
Templates in §4.

### Script layout

```
analysis/<topic>_v<N>/
├── script/
│   ├── 00_config.sh           # paths, env vars, producer toggles (mirrors 0_config.sh idiom)
│   ├── 01_<name>.py           # numbered steps, executable directly
│   ├── 02_<name>.py
│   ├── ...
│   ├── _utils.py              # shared helpers — underscore prefix, never numbered
│   ├── _figure_style.py
│   ├── _profile_plot.py
│   └── run_all.sh             # drives the numbered steps
├── data/                      # ignored — intermediate matrices, h5ads
├── processed/                 # tracked: small text outputs (TSVs, signatures, BEDs)
├── figures/                   # tracked: milestone figures only (project's manifest covers them)
├── qc/                        # tracked: QC tables for the analysis
└── logs/                      # ignored
```

- Numbered steps are **Python by default** for `analysis/` work (vs
  bash for the bulk-seq pipeline), but `00_config.sh` stays bash so it
  can be `source`d by `run_all.sh` and inherited as env vars by step
  scripts.
- `_underscore.py` files are local helpers — never numbered, never executed
  standalone. Promotion has **two tiers**: a helper imported by ≥2 sibling
  analyses **within one project** lifts to the project's `script/local/`
  (shared project-wide but still deposit-specific — e.g. Reed's `_resources.py` /
  `_donor_qc.py` / `_marker_library.py`); only a helper genuinely general
  **across projects** graduates to `scripts/pipeline/tools/<topic>/` (§5).

### Standard Python prelude (numbered steps)

```python
"""NN_<name>.py — <one-line purpose>."""
from __future__ import annotations
import os, sys
from pathlib import Path

PROJ_ROOT = Path(os.environ.get("PROJ_ROOT",
                                 Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(PROJ_ROOT / "script"))
from _utils import ...  # noqa: E402
```

`PROJ_ROOT` env override lets a step be run from anywhere (cron, subagent
worktree, parent driver script) without breaking imports.

### Producer / consumer artifact handoff

When `analysis/A/` produces a table that `analysis/B/` and `analysis/C/`
consume, route through an **env var pinned in B's and C's `00_config.sh`**,
not a hardcoded path:

```bash
# analysis/tf_programs_v2/script/00_config.sh
export FOOTPRINT_TSV="${FOOTPRINT_TSV:-$PROJ_ROOT/../footprint_v2/processed/footprints.tsv}"
```

Lets the consumer flip between producer versions (`footprint/`
vs `_v2/`) by setting `FOOTPRINT_TSV=...` at invocation, without
rewriting paths inside the numbered steps. Pattern in use:
`FOOTPRINT_TSV`, `SKIP_ACTIVE_TF_FILTER`, `PROJ_ROOT`.

---

## §12 — Single-cell / multiome conventions

`seq/multiome/`, `seq/scrna/`, `seq/atac/` (when it's single-cell) use a
different toolchain from the bulk-seq pipelines and need extra discipline
to stay stable on NFS.

### Env

Use the **`sc`** env, not `bio`. `sc` has scanpy, anndata, SnapATAC2,
scrublet, scvi-tools, pyDESeq2, harmonypy; `bio` does not. Snapshots in
`scripts/env/lock/sc.<date>.yml` (same convention as `bio`; §8).

### h5ad / NFS rules

`~/work/` is NFS, and HDF5 file locking on NFS is unreliable:

- Always `export HDF5_USE_FILE_LOCKING=FALSE` before any anndata I/O.
- Stage h5ad > 1 GB to local `/tmp` for writing, then `mv` back.
- **Never write the same NFS h5ad from two processes** — silent
  corruption. Check for orphan subagents before re-running a write step.
- SnapATAC2 per-donor reads: open with `backed='r'` (string, not bool);
  CSR row-slice is ~150× faster than column-slice; barcodes are plain
  (no donor prefix).

### Naming discipline

- Cell-state / subtype labels must not collide with established vocabulary from a
  **different** concept (e.g. immune-cell polarization terms); a collision
  misleads readers. Pick project-specific names that don't overload standard terms.
- Each project's README defines its canonical label space; downstream analyses
  **read** those labels, never re-fit them. The specific per-cohort label
  conventions live with the project (and in the private lab notes), not here.

### Statistical defaults

- **No per-cell Mann-Whitney p-values for inference.** Cells are
  pseudoreplicates of donors; per-cell tests p-hack via N inflation
  (10k cells from 40 donors are not 10k independent observations). Use
  **donor-stratified bootstrap CI on a per-donor summary** (mean log1p
  expression, fraction expressing, AUCell median) instead. Generalizes
  to any pseudoreplicated data, but in practice this is where it bites.
- **Cell-weighted primary, donor-weighted sanity.** Report the
  cell-weighted statistic as the primary result; show the donor-weighted
  version alongside as a sanity check. Discrepancy = composition
  confound; investigate before publishing.
- **Dotplots default to absolute scale.** `sc.pl.dotplot(standard_scale=None, cmap="Reds")`.
  `standard_scale="var"` column-z-scores; broadly-expressed markers
  look "hollow" because 1.1× differences get stretched to full-range contrast.

### Parallel compute (joblib / sklearn / pyDESeq2)

BLAS thread caps set after `import numpy` are no-ops. Two rules:

- Export `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS`
  in the shell *before* launching the Python step (or set them at the
  very top of the script, before any scientific imports).
- Inside joblib workers, use `threadpoolctl.threadpool_limits(N)`
  explicitly — worker processes inherit a fresh BLAS that ignores the
  parent's env vars.

`pyDESeq2 n_cpus=4` whenever ≥2 chains run in parallel (combined with the
16-core system cap; §13).

---

## §13 — System compute cap

linux01 has **16 CPU cores**. Cap *total* concurrent CPU usage across all
parallel jobs (pyDESeq2 `n_cpus`, joblib `n_jobs`, BLAS thread pools,
parallel chains, subagents doing heavy work) at 16. Monitor with
`uptime`; kill stale long-running jobs proactively when load goes above
nominal.
