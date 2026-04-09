# CUT&RUN Upstream Pipeline

Standard pipeline for CUT&RUN data: trimming → alignment → HOMER tags → peak calling → QC.

## Workflow

```
Raw FASTQs  →  1_trim_qc  →  2_bowtie2  →  3_homer_tags  →  4.1_peak_macs3  →  5_qc
   (data/)       (cleandata/)    (align/bam)    (align/tags)   4.2_peak_homer     (multiqc)
                                                              (peaks/macs3, peaks/homer)
```

## Steps

| Step | Script | Description |
|------|--------|-------------|
| 0 | `0_config.sh` | Central config: paths, genome, params. **Edit this first.** |
| 1 | `1_trim_qc.sh` | Trim Galore + FastQC |
| 2 | `2_bowtie2.sh` | Bowtie2 alignment, bigWig tracks, BAM QC |
| 3 | `3_homer_tags.sh` | HOMER tag directories from BAMs |
| 4.1 | `4.1_peak_macs3.sh` | MACS3 peak calling (narrow/broad) + blacklist filter |
| 4.2 | `4.2_peak_homer.sh` | HOMER findPeaks + annotatePeaks + blacklist filter |
| 5 | `5_qc.sh` | preseq, phantompeakqualtools, plotFingerprint, PCA/correlation, MultiQC |

## Usage

```bash
# Run full pipeline
./run_all.sh

# Run only some steps (override via env)
RUN_TRIM=0 RUN_BOWTIE2=0 RUN_PEAK_MACS3=1 RUN_PEAK_HOMER=1 ./run_all.sh

# Run only QC
RUN_TRIM=0 RUN_BOWTIE2=0 RUN_HOMER_TAGS=0 RUN_PEAK_MACS3=0 RUN_PEAK_HOMER=0 RUN_QC=1 ./run_all.sh
```

## Setup

1. **Copy** this folder into your project, e.g. `script/` (same layout as other CUTRUN projects).
2. **FASTQs** — Place paired-end `*_R1.fastq.gz` / `*_R2.fastq.gz` under **`${BASE}/data/`** (default `RAW_DIR`). Symlinks are fine. To download from IGM / ENA or merge Illumina lanes first, use **[`../tools/prep/`](../tools/prep/README.md)** from the repo (`merge_lanes_inplace.sh`, `link_merged_fastqs.sh`, `link_fastq.sh`, etc.).
3. **Edit** `0_config.sh`:
   - `BASE` → project root (e.g. `$HOME/work/seq/CUTRUN/my_project`)
   - `RAW_DIR` → usually `${BASE}/data` (trim input)
   - `CONDA_BIO_ENV` — optional; defaults to `$HOME/miniforge3/envs/bio` and prepends `bin` to `PATH`. Set to `""` to skip, or point at another env.
   - `GENOME` / `GENOME_INDEX` → Bowtie2 index path
   - `BLACKLIST` → blacklist BED path
   - Step toggles **`RUN_*`** — defaults in `0_config.sh`; override per run, e.g. `RUN_TRIM=0 RUN_BOWTIE2=0 RUN_PEAK_MACS3=1 ./run_all.sh` (variables are read from the environment when `0_config.sh` is sourced).
4. **Create** `peakcall_groups.tsv` at project root (path set by `PEAKCALL_GROUPS_FILE`): `ip_bam<TAB>control_bam<TAB>name<TAB>type`

## Requirements

- **Tools** (on `PATH`, or via `CONDA_BIO_ENV`): trim_galore, cutadapt, fastqc, bowtie2, samtools, deeptools (`bamCoverage`), homer, bedtools, macs3 (if `RUN_PEAK_MACS3=1`), preseq, phantompeakqualtools (`run_spp.R`), multiqc
- **R** + packages used by `5_qc.sh` / MultiQC as needed on your system
