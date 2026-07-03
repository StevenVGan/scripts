# __PROJECT__

**Assay:** __ASSAY__
**Genome:** __GENOME__
**SE/PE:** __LAYOUT__
**Owner / submitter:** __OWNER__
**Date received / batch:** __DATE__
**Status:** active
**Upstream source:** __UPSTREAM__

## What this project is
TODO — 2–3 sentences: biological question, IP target(s), the comparison being made.

## Key analyses
- (none yet — add `analysis/<name>_v<N>/` rows as they appear)

## Run
- Pipeline: `scripts/pipeline/__ASSAY_DIR__/` (bulk-seq); env lock recorded in `references.tsv` after a run.
- Standard invocation: `cd script && ./run_all.sh`  — `BASE` auto-derives from the project root, no edit needed.
- Resume from a step: `RUN_TRIM=0 RUN_BOWTIE2=0 ./run_all.sh` (see `0_config.sh` for `RUN_*` toggles).
