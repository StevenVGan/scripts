# Pipeline

Main CUT&RUN analysis pipeline and shared utilities.

## Contents

| Folder | Purpose |
|--------|---------|
| **cutrun/** | Standard CUT&RUN upstream pipeline: trim → align → peak call → QC |
| **tools/** | Reusable utilities: heatmaps, BED liftover, merge peaks, link FASTQs, subsample |

## cutrun

Editable, config-driven pipeline. Copy into each project and adjust `0_config.sh`.

See [cutrun/README.md](cutrun/README.md) for step-by-step details.

## tools

Standalone scripts. Most have config blocks at the top; edit paths before running.

See [tools/README.md](tools/README.md) for usage of each tool.
