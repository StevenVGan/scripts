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

1. **Copy** this folder into your project, e.g. `script/cutrun/` or `script/`
2. **Edit** `0_config.sh`:
   - `BASE` → project root (e.g. `$HOME/work/seq/CUTRUN/my_project`)
   - `RAW_DIR` → where raw/symlinked FASTQs are
   - `GENOME` / `GENOME_INDEX` → Bowtie2 index path
   - `BLACKLIST` → blacklist BED path
3. **Create** `peakcall_groups.tsv` (format: `ip_bam\tcontrol_bam\tname\ttype`)

## Requirements

- trim_galore, cutadapt, fastqc, bowtie2, samtools, deeptools (bamCoverage), homer, bedtools, macs3
- Run inside conda env with these tools on PATH
