# Tools

Reusable utilities for sequencing analysis. **Bash** helpers (`*.sh`) usually use a **config block** at the top—edit paths before running. **Python** tools (`go_enrichr.py`, `annotation_pie.py`) are **CLI-driven** (`--help` for flags); run with `python3 ./script.py` or `./script.py` if executable.

**Working directory:** run **`./scriptname.sh`** / **`python3 ./script.py`** examples from **`scripts/pipeline/tools/`** unless noted (e.g. **`prep/`** examples use `cd prep`).

## Prep (upstream FASTQs)

**[prep/](prep/)** — IGM **`download_fastq`**, ENA **`download_geo_fastq_ena`**, Illumina **`link_fastq`**, **`merge_lanes_inplace`**, **`link_merged_fastqs`**. See **[prep/README.md](prep/README.md)**.

## Topic subfolders

Methodology families with multiple related scripts get their own subfolder
under `tools/`. Convention is documented in
[../../CONVENTIONS.md](../../CONVENTIONS.md) §5. Examples (created when
≥2 related scripts otherwise sit flat at the top of `tools/`):

- `cleavage/` — single-base CUT&RUN cleavage methodology
  (BAM→cut bigwig, motif cutoff sweeps, footprint plots)
- `cobinding/` — ChIP/CUT&RUN cobinding consensus + downstream
- *...add as families emerge*

The family name `<topic>/` matches the prefix of joint repos that consume
it: e.g. `tools/cleavage/` ↔ `seq/_joint/cleavage_sites_*/`. Flat
single-script utilities (`heatmap.sh`, `lift_bed.sh`, `go_enrichr.py`)
stay at the top level.

## ChIP downstream reference — moved

**[chip_downstream_reference/](chip_downstream_reference/)** is now a
**one-page README pointer** to `seq/_joint/MCF7_ER_p65_cobinding/`, the
joint analysis repo that owns the (formerly duplicated) GO Enrichr,
annotation composition, and signal-profile scripts. Single source of
truth — no more edit-here-copy-there sync. See
[chip_downstream_reference/README.md](chip_downstream_reference/README.md).

For generic bigWig heatmaps (peaks / TSS / gene bodies) without Method 2 site-class logic, use **heatmap.sh** in the table below.

**Generic downstream CLIs** (Enrichr + annotation pies; CLI-style like `peak_ops.sh`):

- **[go_enrichr.py](go_enrichr.py)** — gene list, HOMER annotate table (`Gene Name`), or BED name column → Enrichr TSV + top-terms bar chart (`--gene-set` e.g. `GO_Biological_Process_2023`).
- **[annotation_pie.py](annotation_pie.py)** — HOMER `Annotation` column (or arbitrary TSV column / one string per line) → single pie chart PNG.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| **heatmap.sh** | deepTools heatmaps of bigWig tracks over peaks/TSS/genes | Edit config (BASE, BED_FILES, BW_FILES, REGION), then `./heatmap.sh` |
| **go_enrichr.py** | Enrichr enrichment + bar plot from genes / annotatePeaks / BED names | `./go_enrichr.py --input FILE --source genes\|annotate\|bed --out-prefix OUT/prefix [--gene-set NAME] [--gene-column COL]` |
| **annotation_pie.py** | Pie chart of genomic annotation categories | `./annotation_pie.py --input FILE [--format auto\|homer\|tsv\|lines] -o plot.png [--annotation-column COL]` |
| **lift_bed.sh** | LiftOver BED files between assemblies (e.g. hg19→hg38) | `./lift_bed.sh INPUT FROM TO OUTPUT [CHAIN]` |
| **merge_peaks.sh** | Merge multiple BED/HOMER peak files into one | `./merge_peaks.sh OUTPUT.bed INPUT1 [INPUT2 ...]` |
| **intersect_peaks.sh** | Wrapper for `peak_ops.sh --mode intersect` (same BED + UpSet/Venn) | `./intersect_peaks.sh [--viz …] [--slop BP] [--names "N1,N2,..."] OUTPUT.bed INPUT1 INPUT2 [...]` |
| **peak_ops.sh** | Peak set ops: intersect, distinct, or union + optional UpSet / Venn | `./peak_ops.sh --mode MODE [--viz none\|upset\|venn\|both] [--slop BP] [--names "N1,N2,..."] OUTPUT.bed INPUT1 INPUT2 [...]` |
| **getfasta.sh** | Extract sequences from reference genome for BED regions | `./getfasta.sh [--slop BP] [--no-center] [--tab] OUTPUT.txt INPUT1.bed [INPUT2.bed ...]` |
| **subsample_data.sh** | Subsample FASTQs (e.g. 1M reads) for testing | Uses `0_config.sh` from `../cutrun/` or set `CONFIG_FILE=/path/to/project/0_config.sh`; needs `seqtk` |

## Examples

**Heatmap**
```bash
# Edit BASE, BED_FILES, BW_FILES, REGION (Peaks|TSS|Genes), OUTPUT_NAME
./heatmap.sh
```

**Enrichr GO (generic)**

```bash
# One gene symbol per line
./go_enrichr.py --input genes.txt --source genes --out-prefix results/my_peaks

# HOMER annotatePeaks.txt (uses Gene Name column by default)
./go_enrichr.py -i peaks.annotatePeaks.txt --source annotate --out-prefix results/peaks \
  --gene-set GO_Biological_Process_2023 --organism human

# BED: column 4 = gene symbol per interval
./go_enrichr.py -i regions.bed --source bed --bed-gene-col 4 --out-prefix results/bedgenes --top-n 12
```

Requires: `gseapy`, `matplotlib`, network access to Enrichr. TSV: `*_enrichr.tsv`; plot: `*_top_terms.png` (omit with `--no-plot`). Plotting picks **Term** and an adjusted/nominal **P-value** column case-insensitively when possible (`Adjusted P-value` preferred).

**Annotation composition pie**

```bash
# HOMER annotatePeaks (auto-detects Annotation column)
./annotation_pie.py -i rep1.annotatePeaks.txt --format homer -o annotation.png --title "My peaks"

# Raw annotation strings, one per line
./annotation_pie.py -i ann_strings.txt --format lines -o pie.png --counts-tsv counts.tsv

# Same as --format homer when the first line is a TSV header containing Annotation
./annotation_pie.py -i rep1.annotatePeaks.txt --format auto -o annotation_auto.png
```

**Lift BED**
```bash
./lift_bed.sh ./MCF7_Amir_hg19 hg19 hg38 ./MCF7_Amir_hg38
./lift_bed.sh my.bed hg19 hg38 ./out
```

**Merge peaks**
```bash
./merge_peaks.sh merged.bed rep1_peaks.bed rep2.annotatePeaks.txt
```

**Intersect peaks** (`intersect_peaks.sh` → `peak_ops.sh --mode intersect`; optional `--viz` like peak_ops)
```bash
./intersect_peaks.sh out.bed rep1.annotatePeaks.txt rep2.annotatePeaks.txt
./intersect_peaks.sh --viz both --slop 250 --names "R1,R2,R3,R4" peaks/ERa_intersect.bed peaks/*.bed
```

**Peak set operations** (intersect / distinct / union; plots controlled by `--viz`)

- **`--viz upset`** (default): UpSet PDF next to output (`*_upset.pdf`), same logic as [peaks_ops_upsetR.R](peaks_ops_upsetR.R) (partition counts for set sizes).
- **`--viz venn`**: VennDiagram PDF (`*_venn.pdf`) for **2–4** input peak sets only; skipped with a message when **n > 4**.
- **`--viz both`**: UpSet and Venn when n ≤ 4; UpSet only when n > 4.
- **`--viz none`**: BED only, no R.
- **`intersect_peaks.sh`** is a thin `exec` wrapper around `peak_ops.sh --mode intersect` (same BED, UpSet, and Venn behavior).

```bash
# intersect: regions in ALL inputs; UpSet + quad Venn (4 samples)
./peak_ops.sh --mode intersect --viz both --slop 250 --names "Veh,E2,TNF,E2TNF" peaks/panel.bed peaks/*.bed

# intersect, plots off (faster batch)
./peak_ops.sh --mode intersect --viz none out.bed rep1.bed rep2.bed

# distinct: regions in exactly the specified sets (partition: in A&B, NOT in C or D)
./peak_ops.sh --mode distinct --viz upset --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/distinct.bed peaks/*.annotatePeaks.txt

# union: regions in ANY input
./peak_ops.sh --mode union --viz upset --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/union.bed peaks/*.annotatePeaks.txt
```
Requires: `bedtools`; R + **UpSetR** and (for Venn) **VennDiagram** (`install.packages(c('UpSetR','VennDiagram'), repos='https://cloud.r-project.org')`). Ensure `Rscript` is on `PATH` (e.g. conda env).

**Get sequences from BED regions**
```bash
# Default: center regions, extend ±500 bp (1001 bp total), output sequences only
./getfasta.sh seqs.txt peaks.bed

# Use original coordinates, no slop
./getfasta.sh --no-center --slop 0 seqs.txt peaks.bed

# Tab format (name, seq) for downstream use
./getfasta.sh --tab seqs.txt rep1.bed rep2.bed
```

**Link FASTQs** (`prep/link_fastq.sh`)

- **Map file (`MAP_FILE`):** each non-comment line: `prefix` TAB `newname` (optional extra columns ignored). Example: Illumina lane prefix or SRR accession → pipeline sample base name (`sample_R1.fastq.gz` / `sample_R2.fastq.gz` in `DEST_DIR`).
- **Behavior:** `ln -sfn` for atomic refresh; `bash` `nullglob` when expanding `prefix*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`; if PE R1 pattern matches nothing, falls back to `${RAW_DIR}/${prefix}.fastq.gz` for SE GEO/ENA downloads. If R2 is missing, removes `${newname}_R2.fastq.gz` in `DEST_DIR` so old symlinks do not linger.
- **Log:** `${MAP_FILE%.tsv}.log.tsv` (e.g. `link_sample.log.tsv` next to `link_sample.tsv`).

```bash
# Option A: edit defaults at top of prep/link_fastq.sh
./prep/link_fastq.sh

# Option B: project wrapper (recommended): export and exec
export RAW_DIR=~/work/raw_seq/my_run/fastq
export DEST_DIR=~/work/seq/my_project/data
export MAP_FILE=~/work/seq/my_project/link_sample.tsv
exec bash ~/work/scripts/pipeline/tools/prep/link_fastq.sh
```

**Multi-lane IGM → merged names → `data/`:** see [prep/README.md](prep/README.md) (`merge_lanes_inplace.sh`, `link_merged_fastqs.sh`).

**ENA / GEO-style SRR downloads** (`prep/download_geo_fastq_ena.sh`; campus IGM FTP is **`prep/download_fastq.sh`**)
```bash
cd prep
# Edit DEST_DIR, SRR_LIST_FILE in download_geo_fastq_ena.sh (replace MY_GEO_RUN), or:
DEST_DIR=~/work/raw_seq/my_study/fastq \
SRR_LIST_FILE=~/work/raw_seq/my_study/srr_accessions.txt \
DOWNLOAD_JOBS=8 \
  ./download_geo_fastq_ena.sh

# Optional: pass list path as first argument (DEST_DIR still from config/env)
./download_geo_fastq_ena.sh /path/to/srr_accessions.txt
```
Requires: `curl`, `md5sum`. Writes `*.part.*` then renames to `SRR....fastq.gz`. Logs and `ena_manifest.tsv` go under `LOG_DIR` (default: parent of `DEST_DIR` + `/logs`). Checksums: `MD5_FILE` (default: parent of `DEST_DIR` + `/md5sum_ena.txt`). **SE libraries only**; paired-end runs with two URLs in `fastq_ftp` are rejected until extended.

## Notes

- **chip_downstream_reference**: now a pointer to
  `seq/_joint/MCF7_ER_p65_cobinding/` (the joint repo is the single source
  of truth for cobinding scripts). Generic **go_enrichr.py** /
  **annotation_pie.py** above are independent — use them directly.
- **prep/** (**download_fastq**, **download_geo_fastq_ena**, **link_fastq**, **merge_lanes_inplace**, **link_merged_fastqs**) and **subsample_data**: set paths via script config or env; project wrappers should `exec` scripts under **`tools/prep/`** where appropriate.
- **peak_ops `--viz`**: Venn diagrams use the same **mutually exclusive peak partitions** as the UpSet right-bar counts ([peak_vennDiagram.R](peak_vennDiagram.R)); not supported for more than four sets (use UpSet).
- **UpSet (5+ sets):** [peaks_ops_upsetR.R](peaks_ops_upsetR.R) omits per-bar colored **queries** when **n > 4** to avoid an UpSetR bug; counts and default UpSet bars are unchanged.
- **subsample_data**: Uses `../cutrun/0_config.sh` when run from `pipeline/tools/`, or set `CONFIG_FILE=/path/to/project/script/0_config.sh` for your project.
