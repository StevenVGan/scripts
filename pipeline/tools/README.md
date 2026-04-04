# Tools

Reusable utilities for sequencing analysis. Each script has a config block at the top; edit paths before running.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| **heatmap.sh** | deepTools heatmaps of bigWig tracks over peaks/TSS/genes | Edit config (BASE, BED_FILES, BW_FILES, REGION), then `./heatmap.sh` |
| **lift_bed.sh** | LiftOver BED files between assemblies (e.g. hg19→hg38) | `./lift_bed.sh INPUT FROM TO OUTPUT [CHAIN]` |
| **merge_peaks.sh** | Merge multiple BED/HOMER peak files into one | `./merge_peaks.sh OUTPUT.bed INPUT1 [INPUT2 ...]` |
| **intersect_peaks.sh** | Intersect peaks (regions in ALL inputs) + UpSet plot | `./intersect_peaks.sh [--slop BP] [--names "N1,N2,..."] OUTPUT.bed INPUT1 INPUT2 [...]` |
| **peak_ops.sh** | Peak set ops: intersect, distinct, or union + UpSet plot | `./peak_ops.sh --mode MODE [--slop BP] [--names "N1,N2,..."] OUTPUT.bed INPUT1 INPUT2 [...]` |
| **getfasta.sh** | Extract sequences from reference genome for BED regions | `./getfasta.sh [--slop BP] [--no-center] [--tab] OUTPUT.txt INPUT1.bed [INPUT2.bed ...]` |
| **link_fastq.sh** | Symlink raw FASTQs into project `data/` from a TSV map (Illumina `*_R1/_R2_001` or SE `SRR.fastq.gz`) | Edit config or `export RAW_DIR DEST_DIR MAP_FILE`, then `./link_fastq.sh` |
| **download_geo_fastq_ena.sh** | Download SRR FASTQs from ENA (HTTPS/curl), parallel jobs | Edit `DEST_DIR`, `SRR_LIST_FILE` (or env overrides); SE-only per SRR |
| **subsample_data.sh** | Subsample FASTQs (e.g. 1M reads) for testing | Uses `0_config.sh` from `../cutrun/` or set `CONFIG_FILE=/path/to/project/0_config.sh`; needs `seqtk` |

## Examples

**Heatmap**
```bash
# Edit BASE, BED_FILES, BW_FILES, REGION (Peaks|TSS|Genes), OUTPUT_NAME
./heatmap.sh
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

**Intersect peaks** (regions in ALL inputs; generates UpSet PDF)
```bash
./intersect_peaks.sh out.bed rep1.annotatePeaks.txt rep2.annotatePeaks.txt
./intersect_peaks.sh --slop 250 --names "R1,R2,R3,R4" peaks/ERa_intersect.bed peaks/*.annotatePeaks.txt
```

**Peak set operations** (intersect / distinct / union; generates UpSet PDF)
```bash
# intersect: regions in ALL inputs
./peak_ops.sh --mode intersect --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/intersect.bed peaks/*.annotatePeaks.txt

# distinct: regions in exactly the specified sets (partition: in A&B, NOT in C or D)
./peak_ops.sh --mode distinct --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/distinct.bed peaks/*.annotatePeaks.txt

# union: regions in ANY input
./peak_ops.sh --mode union --slop 250 --names "ICI_rep1,ICI_rep2,E2_rep1,E2_rep2" peaks/union.bed peaks/*.annotatePeaks.txt
```
Requires: R + UpSetR (`install.packages('UpSetR', repos='https://cloud.r-project.org')`)

**Get sequences from BED regions**
```bash
# Default: center regions, extend ±500 bp (1001 bp total), output sequences only
./getfasta.sh seqs.txt peaks.bed

# Use original coordinates, no slop
./getfasta.sh --no-center --slop 0 seqs.txt peaks.bed

# Tab format (name, seq) for downstream use
./getfasta.sh --tab seqs.txt rep1.bed rep2.bed
```

**Link FASTQs** (`link_fastq.sh`)

- **Map file (`MAP_FILE`):** each non-comment line: `prefix` TAB `newname` (optional extra columns ignored). Example: Illumina lane prefix or SRR accession → pipeline sample base name (`sample_R1.fastq.gz` / `sample_R2.fastq.gz` in `DEST_DIR`).
- **Behavior:** `ln -sfn` for atomic refresh; `bash` `nullglob` when expanding `prefix*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`; if PE R1 pattern matches nothing, falls back to `${RAW_DIR}/${prefix}.fastq.gz` for SE GEO/ENA downloads. If R2 is missing, removes `${newname}_R2.fastq.gz` in `DEST_DIR` so old symlinks do not linger.
- **Log:** `${MAP_FILE%.tsv}.log.tsv` (e.g. `link_sample.log.tsv` next to `link_sample.tsv`).

```bash
# Option A: edit defaults at top of link_fastq.sh
./link_fastq.sh

# Option B: project wrapper (recommended): export and exec unchanged tool
export RAW_DIR=~/work/raw_seq/my_run/fastq
export DEST_DIR=~/work/seq/my_project/data
export MAP_FILE=~/work/seq/my_project/link_sample.tsv
exec bash ~/work/scripts/pipeline/tools/link_fastq.sh
```

**ENA / GEO-style SRR downloads** (campus IGM FTP is `download_fastq.sh` instead)
```bash
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

- **link_fastq**, **download_geo_fastq_ena**, and **subsample_data**: project-specific; edit config paths or set env vars before running. Prefer thin `script/link_fastq.sh` wrappers that export `RAW_DIR` / `DEST_DIR` / `MAP_FILE` and `exec` the tool script.
- **subsample_data**: Uses `../cutrun/0_config.sh` when run from `pipeline/tools/`, or set `CONFIG_FILE=/path/to/project/script/0_config.sh` for your project.
