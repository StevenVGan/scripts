# scripts/maintenance — inactive-project archival

One-time (and re-runnable) cleanup that recovers disk space from sequencing
projects that aren't being actively analyzed, without losing anything that
can't be regenerated from `~/work/raw_seq/` and the pipeline.

## What this does, per project

For every row in `projects.tsv` whose `status` column is `cleanup`, the
script `archive_inactive.sh` performs three actions:

1. **BAM → CRAM (lossless).** Each `*.bam` under `align/bam/` (or `align/`
   for legacy projects) is converted to CRAM against the project's reference
   genome, indexed, and verified (`samtools view -c` parity + `samtools
   quickcheck`). Only after both checks pass is the original BAM (and `.bai`)
   removed.
2. **Delete `align/tags/`.** HOMER tag directories are regeneratable.
3. **Delete `cleandata/` and `cleandata2/`.** Trimmed FASTQs are
   regeneratable.

Bigwigs (`align/track/*.bw`), peaks, multiqc, raw symlinks, scripts, logs,
and the entire `align/bamqc/` tree are **never** touched. Active projects
and `keep`-status projects are skipped entirely.

## What changes vs what's preserved

| Path | Cleanup-status projects | Active / keep projects |
|---|---|---|
| `align/bam/*.bam` (and legacy `align/*.bam`) | **BAM → CRAM** (lossless), original BAM removed after verification | unchanged |
| `align/bam/*.bai` | **removed** with the BAM it indexed | unchanged |
| `align/tags/` | **removed entirely** | unchanged |
| `cleandata/`, `cleandata2/` | **removed entirely** | unchanged |
| `align/track/*.bw` (bigwigs) | **kept** | unchanged |
| `align/bamqc/` | kept | unchanged |
| `peaks/`, `multiqc/`, `data/`, `script/`, `logs/`, `matrix/`, `ori_data/` | kept | unchanged |
| `seq/_joint/*` | n/a (joint analyses are status=active) | unchanged |
| `seq/raw_seq/`, `seq/ref/`, `scripts/`, `~/miniforge3` | out of scope | unchanged |

## Project-status table

`projects.tsv` is the single source of truth. Statuses:

| Status | Meaning |
|---|---|
| `active` | Currently being analyzed (or about to be). **Untouched.** |
| `keep`   | Kept for record/storage reasons (e.g. csRNA single-date stubs, the chip `_input` control ref). **Untouched.** |
| `cleanup` | The script processes this project end-to-end. |

## Running it

The script is config-driven (env vars at the top); no CLI flags. Activate
the `bio` conda env first so `samtools` is on `PATH`. From any directory:

```bash
cd ~/work/scripts/maintenance

# Phase 1 — dry-run inventory (no deletions). Writes inventory.tsv with
# per-project bytes and detected genome.
./archive_inactive.sh inventory

# Phase 2 — smoke test on one mid-size project. Suggested target:
#   seq/cnr/CnR_251204_CTCF_HEK_Priyanka  (~19G, hg38, modern align/bam/ layout)
# Then sanity-check log under logs/, and verify CRAM round-trips back to BAM:
./archive_inactive.sh run seq/cnr/CnR_251204_CTCF_HEK_Priyanka

# Phase 3 — full batch. Loops through every cleanup row in projects.tsv.
./archive_inactive.sh run

# Phase 4 — accounting.
./archive_inactive.sh audit       # writes audit.tsv with before/after/saved
```

Per-project logs land under `scripts/maintenance/logs/`.

### Useful env overrides

```bash
REF_HG38=/some/other/hg38.fa ./archive_inactive.sh run ...
THREADS=16 ./archive_inactive.sh run
DEFAULT_GENOME=hg38 ./archive_inactive.sh run    # used only when auto-detect fails
```

Reference defaults:

- `REF_HG38=/mnt/share/archive/bkup/ref/genome/hg38/hg38.fa`
- `REF_HG19=/mnt/share/archive/bkup/ref/genome/hg19/hg19.fa`
- `REF_MM10=/mnt/share/archive/bkup/ref/genome/mm10/mm10.fa`

### Genome detection

Per project, in order:

1. Parse `GENOME=` from `script/0_config.sh`.
2. If absent, read `chr1` length from the BAM header (`@SQ` line) and match
   to a known assembly — hg38=248956422, hg19=249250621, mm10=195471971.
3. Fall back to `DEFAULT_GENOME` env var.

If none of those work, BAM→CRAM is skipped for that project (with a clear
`WARN` line in the log) and tags + cleandata are still removed.

## Regenerating what was removed

All three deletions are reversible from the standard pipeline plus
`~/work/raw_seq/`. The exact commands depend on the assay; CUT&RUN examples
below — substitute `pipeline/csrna/` or `pipeline/pro/` as needed.

### CRAM → BAM

```bash
# Single file:
samtools view -@ 8 -b -T /mnt/share/archive/bkup/ref/genome/hg38/hg38.fa \
  -o sample_sorted.bam sample_sorted.cram
samtools index sample_sorted.bam

# Whole project, all CRAMs at once:
cd ~/work/seq/<assay>/<project>/align/bam
for c in *.cram; do
  samtools view -@ 8 -b -T <ref.fa> -o "${c%.cram}.bam" "$c"
  samtools index "${c%.cram}.bam"
done
```

CRAM is what MACS3 and HOMER `makeTagDirectory` cannot ingest directly, so
this conversion is the prerequisite for re-running peak calling or
re-tagging on an archived project.

### Tag directory regeneration

From the per-sample BAM (after CRAM→BAM if the project was archived):

```bash
makeTagDirectory align/tags/<sample> align/bam/<sample>_sorted.bam
```

Or sweep all samples by following the project's
`script/3_homer_tags.sh` step (since pipeline scripts are config-driven, the
right invocation is just to re-run that step in the project's `script/`
folder with `RUN_HOMER_TAGS=1` and the other `RUN_*` toggles off).

### Trimmed FASTQ regeneration

The trim step reads from `data/` (raw symlinks) and writes to `cleandata/`.
Re-run by toggling only `RUN_TRIM`:

```bash
cd ~/work/seq/<assay>/<project>/script
RUN_TRIM=1 RUN_BOWTIE2=0 RUN_HOMER_TAGS=0 RUN_PEAK_MACS3=0 \
  RUN_PEAK_HOMER=0 RUN_QC=0 ./run_all.sh
```

For PRO-seq projects the second-pass `cleandata2/` is regenerated by the
same step in `pipeline/pro/`.

## Safety properties

- **Reference is the only external dependency for CRAM.** As long as the
  same reference FASTA is preserved, CRAM is bit-exact for the alignments;
  there is no information loss vs the source BAM. If a verification check
  fails, the script bails on that BAM and leaves it intact.
- **Nothing under `~/work/raw_seq/` or `data/` is touched** — the raw
  FASTQs are the ultimate ground truth for regenerating everything else.
- **Per-project failures are non-fatal** in batch mode; the loop continues
  and reports the failure count at the end.
- **Re-running is safe.** If a CRAM already exists for a given BAM, the
  script skips that pair without re-encoding; tags and cleandata removal
  are idempotent.

## Tuning the cleanup set

Edit `projects.tsv`. Move a path's `status` to:

- `active` — currently being analyzed; never touched.
- `keep` — preserved for record/storage; never touched.
- `cleanup` — processed by the script.

Comment lines aren't supported (the parser is `tail -n +2 | while read`).
The header row is hard-coded (`status\tpath\tnotes`).

If you want to additionally drop `matrix/` (deepTools heatmap matrices —
fully regeneratable from bigwig + BED) or `align/bamqc/` for cleanup
projects, add the rule to `archive_one_body` in `archive_inactive.sh`.
Both targets are intentionally left alone by the current script because the
savings are small relative to the cost of regenerating heatmaps.
