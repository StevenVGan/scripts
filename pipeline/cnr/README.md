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
| 5 | `5_qc.sh` | preseq, plotFingerprint, bamPEFragmentSize, estimateReadFiltering, PCA/correlation, optional phantompeakqualtools (`RUN_SPP=1`), MultiQC |

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

## Parameters worth knowing (set in `0_config.sh`; changed 2026-10-20 after a full parameter review on a PE101 CUT&RUN batch)

| Knob | Default | What it does |
|------|---------|--------------|
| `BT2_PAIR_FLAGS` | `--no-mixed --no-discordant --dovetail -I 10 -X 700` | Concordant-pairs-only alignment (Henikoff / EpiCypher / nf-core recipe). MACS3 BAMPE always used proper pairs only; this makes HOMER tags, bigWigs and FRiP see the same fragments, keeps adapter-remnant dovetails and 500–700 bp di-nucleosome pairs. `BT2_PAIR_FLAGS=""` = the permissive pre-2026-10 behaviour (every project before then). |
| `TRIM_2COLOUR` | `1` | `--2colour $TRIM_QUAL` (poly-G aware, NovaSeq/NextSeq) instead of `-q`. `0` for 4-colour instruments. |
| `TRIM_CPU` / `FASTQC_THREADS` | `4` / `8` | Trim Galore `--cores 4` ≈ 15 threads for PE; 8 was ~27, over the node's 16-CPU cap. |
| `TMPDIR` (override `PIPELINE_TMPDIR`) | `$BASE/scratch/tmp` | Scratch for deepTools bedGraph chunks and MACS3 `--tempdir`, kept out of /tmp: the nightly /tmp cleaner (`guard/clean_tmp.sh`, 04:00) deleted `_deeptools_*.bg` under a running bamCoverage on 2026-10-21. The guard now skips them while deepTools runs; pipeline scratch stays off /tmp regardless. |
| `SORT_CPU` / `SORT_MEM` | `4` / `2G` | `samtools sort -@` counts *additional* threads; 12 + 4 keeps the bowtie2 pipe at the cap. |
| `BAMCOV_IGNOREDUP`, `HOMER_TBP`, `MACS3_KEEPDUP` | `1`, `1`, `1` | Positional-duplicate caps (historical). For CUT&RUN TF arms run a `MACS3_KEEPDUP=all` variant into `PEAK_DIR=$BASE/peaks_keepdupall` and report both. HOMER caps per 5′ end with 0.5-weight mates — see the comment block. |
| `FINDPEAKS_FACTOR_EXTRA` | empty | Extra `findPeaks` args for `-style factor` rows only, e.g. `-size 200` (fixed TF window; HOMER's `auto` estimate is 140–160 bp in working ERα arms but 420–600 bp in low-enrichment arms). |
| `HOMER_CHECKGC` | `1` | `makeTagDirectory -genome -checkGC`: tag-GC vs genome-GC per library, drawn by MultiQC (post-alignment contamination readout). |
| `RUN_SPP` | `0` | phantompeakqualtools off: on paired-end BAMs with both mates NSC/RSC are tautological and it costs hours. |
| `BAMQC_REGION`, `BAMQC_BINSIZE`, `BAMQC_DISTANCE` | empty, `10000`, `90000` | PCA/correlation on ~31k genome-wide spaced bins with blacklist, duplicate and MAPQ ≥ 10 filters (was chr19 only). |
| `TAG_DIR`, `PEAK_DIR`, `BAM_DIR` | derived | Env-overridable so variant calls never overwrite the canonical set (the reuse checks are existence-based). Run variants with `RUN_QC=0`: `MULTIQC_OUT` is not variant-aware, and MultiQC scanning `$BASE` would merge the variant's peak stats into the canonical report under the same sample names. |

Also since 2026-10-20: `preseq -P` on PE BAMs; `bamPEFragmentSize` + `estimateReadFiltering` tables; per-control `plotFingerprint --JSDsample` metrics for every distinct control in `peakcall_groups.tsv` (tables in `align/bam/bamqc/fingerprint_metrics_<control>.txt`, kept out of MultiQC because sample names collide across controls); FRiP on the proper-pair universe (`FRiP_properpair`), IP BAM looked up by row name so contrast / `_noctrl` / `_narrow` rows resolve; MACS3 outputs of a row are removed before each call; HOMER `findPeaks` tables are gzipped instead of deleted; MultiQC shows properly-paired %, MQ0 %, insert size and the MACS3 redundancy columns, and ignores `qc/spikein/` and `qc/contam/`. For a sharp histone mark (H3K4me3) add a duplicate row typed `TF` (`<name>_narrow`) to get a MACS3 narrowPeak set next to the `Histone` (broad / HOMER histone-style) call.

**Forks:** `../atac/`, `../csrna/`, `../pro/` carry their own copies of these steps. They received only the `TMPDIR` scratch handling (2026-10-21); the other 2026-10-20 parameter decisions were NOT ported.

## Requirements

- **Tools** (on `PATH`, or via `CONDA_BIO_ENV`): trim_galore, cutadapt, fastqc, bowtie2, samtools, deeptools (`bamCoverage`, `plotFingerprint`, `multiBamSummary`, `bamPEFragmentSize`, `estimateReadFiltering`), homer, bedtools, macs3 (if `RUN_PEAK_MACS3=1`), preseq, multiqc; phantompeakqualtools (`run_spp.R`) only if `RUN_SPP=1`
- **R** + packages used by `5_qc.sh` / MultiQC as needed on your system
