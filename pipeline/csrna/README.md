# csRNA / nascent RNA upstream pipeline

Fork of [../cutrun/README.md](../cutrun/README.md) for **stranded** nascent/csRNA-style libraries: optional **second-pass cutadapt** on Trim Galore outputs, optional **MultiQC immediately after trimming**, extra **Trim Galore** / **Bowtie2** args from config, **strand-separated bigWigs**, and HOMER **`makeTagDirectory -sspe`** when enabled.

## Workflow

```
Raw FASTQs (data/)  →  1_trim_qc  →  1c_trim_pass2 (optional)  →  1b_trim_multiqc (optional)  →  2_bowtie2  →  3_homer_tags  →  4.1 / 4.2  →  5_qc
```

## Default project

`0_config.sh` defaults `BASE` to the combined HepG2 project:

`$HOME/work/seq/PROseq/240216_240701_csRNA_HepG2-ER3XHA_Steven`

Override: `BASE=/path/to/project ./run_all.sh`

## Key variables (`0_config.sh` or env)

| Variable | Default | Meaning |
|----------|---------|---------|
| `RUN_POST_TRIM_MQC` | `1` | Run `1b_trim_multiqc.sh` after trimming; set `0` once adapters look clean |
| `RUN_TRIM_PASS2` | `1` | Run `1c_trim_pass2.sh` on `*_val_*.fq.gz` / `*_trimmed.fq.gz`, then refresh trim FastQC (set `0` to skip after QC is stable) |
| `PASS2_MODE` | `trim_galore` | `trim_galore` = second Trim Galore (Illumina-oriented adapter autodetection); `cutadapt` = explicit `PASS2_ADAPTER_3P` / optional 5′ (most reproducible) |
| `TRIM_GALORE_PASS2_EXTRA` | empty | Extra Trim Galore args for pass-2 only (e.g. `--illumina`), space-separated |
| `PASS2_ADAPTER_3P` | TruSeq universal | Used when `PASS2_MODE=cutadapt`: 3′ adapter for cutadapt `-a`/`-A` (PE) |
| `PASS2_ADAPTER_5P` | empty | Optional 5′ adapter (`-g` / `-G`); set from FastQC overrepresented sequences if needed |
| `PASS2_ADAPTER_5P_R2` | empty | R2-only 5′ sequence when it differs from `PASS2_ADAPTER_5P` |
| `PASS2_MINLEN` | same as `TRIM_MIN_LENGTH` | Minimum read length after pass-2 |
| `PASS2_QUAL` | same as `TRIM_QUAL` | Quality trim threshold in pass-2 |
| `PASS2_OVERLAP` | `4` | cutadapt minimum overlap |
| `PASS2_ERROR_RATE` | `0.15` | cutadapt max error rate |
| `TRIM_GALORE_EXTRA` | empty | Space-separated extra args to Trim Galore (tune after raw FastQC) |
| `BT2_EXTRA` | empty | e.g. `--fr` or `--rf` for stranded PE orientation |
| `STRAND_BIGWIG` | `1` | Emit `_fwd.bw` and `_rev.bw` (`bamCoverage --filterRNAstrand`) |
| `COMBINED_BIGWIG` | `1` | When `STRAND_BIGWIG=1`, also emit unstranded `${sample}.bw` |
| `HOMER_SS_PE` | `1` | Pass `-sspe` to `makeTagDirectory` (disable if chemistry differs) |
| `TRIM_MIN_LENGTH` | `20` | Post-trim minimum length (HOMER csRNA example uses ~20 nt floor) |
| `TRIM_QUAL` | `25` | Phred threshold for 3′ quality trimming (20 is a common looser baseline) |
| `BAMCOV_IGNORE_DUP` | `0` | Set `1` to add `--ignoreDuplicates` to `bamCoverage` (HOMER discourages dedup for csRNA signal) |

## Parameter notes (literature / docs)

- **Trimming / length:** [HOMER csRNA-seq tutorial](http://homer.ucsd.edu/homer/ngs/csRNAseq/index.html) states csRNA captures **~20–65 nt** capped RNAs and recommends trimming 3′ adapters; its `homerTools trim` example discards reads **shorter than 20 nt** after trim. Reads **under ~15–20 nt** are described as largely unusable for alignment. [Trim Galore](https://github.com/FelixKrueger/TrimGalore) wraps Cutadapt; for **Illumina smallRNA** kits it can lower the default minimum length automatically—generic PE defaults are higher, so an explicit **`--length` / `TRIM_MIN_LENGTH`** matters for csRNA.
- **Pass-2:** Short-insert PE can leave **residual adapter** (often visible on R2 in FastQC). With **`RUN_TRIM_PASS2=1`** (default), **`PASS2_MODE=trim_galore`** runs a **second Trim Galore** pass (same adapter heuristics as step 1) on the `*_val_*` reads, then rebuilds `cleandata/fastqc/` so **`1b`** MultiQC reflects the final state. Pass-2 Trim Galore logs and reports are kept as `*_pass2_*` in `TRIM_DIR`. For maximum control or nonstandard kits, set **`PASS2_MODE=cutadapt`** and tune **`PASS2_ADAPTER_3P`**. Set **`RUN_TRIM_PASS2=0`** after QC is stable if you want to skip pass-2 on reruns.
- **Quality:** Phred **≥20** is a common minimum-quality convention; **`TRIM_QUAL=25`** is stricter (fewer bases trimmed at 3′ end if quality is good).
- **Mapping:** HOMER suggests a **splicing-aware** aligner (**STAR** / **HISAT2**) for csRNA in case of spliced contaminant RNAs; this fork still uses **Bowtie2** (CUT&RUN-like). **`--very-sensitive`** is a reasonable general preset; **`--very-sensitive-local`** is sometimes preferred when reads still carry low-quality or soft ends (discussed on forums for difficult reads—test both if alignment rate is low). Set **`BT2_EXTRA`** for stranded orientation (**`--fr`** / **`--rf`**).
- **Coverage / normalization:** [deepTools `bamCoverage`](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html): **`RPGC`** + **`--effectiveGenomeSize`** gives 1× genome scaling (appropriate with hg38 size table). **Do not use `--extendReads`** for spliced/RNA-style mapping (default off here). **`--filterRNAstrand`** assumes **TruSeq/dUTP-style** stranded libraries (read2 in RNA direction); flip interpretation if your prep is the opposite ([deepTools note](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html)).
- **Duplicates:** HOMER’s csRNA page recommends **not** removing PCR duplicates for routine csRNA quantification (dynamic range); **`BAMCOV_IGNORE_DUP=0`** matches that for bigWigs. **`5_qc.sh`** still uses `--ignoreDuplicates` for `plotFingerprint` (fingerprint-style QC only).
- **HOMER tags:** Tutorial mentions **`-maxlen 65`** etc. to filter contaminants; not wired by default—extend `makeTagDirectory` manually if needed. Optional **STAR** alignment is a larger pipeline change.

## `peakcall_groups.tsv`

Same four-column tab format as CUT&RUN (see `seq/CUTRUN/260401_CnR_ERa_OGG1_MCF7_Priyanka/peakcall_groups.tsv`). Place at `${BASE}/peakcall_groups.tsv`. The combined HepG2 project ships an example at its root.

## Running under `screen`

```bash
screen -S csRNA_hepg2
cd ~/work/seq/PROseq/240216_240701_csRNA_HepG2-ER3XHA_Steven/scripts/csRNA
# Pass-2 is on by default; set stranded PE orientation from your prep, e.g.:
BT2_EXTRA="--fr" ./run_all.sh
```

**Trim and QC only** (no alignment or downstream):  
`RUN_TRIM=1 RUN_BOWTIE2=0 RUN_HOMER_TAGS=0 RUN_PEAK_HOMER=0 RUN_PEAK_MACS3=0 RUN_QC=0 ./run_all.sh`  
(add `RUN_TRIM_PASS2=0` if you want trim-only without pass-2)  
(still runs **`1b`** if `RUN_POST_TRIM_MQC=1`).

After validating trimming, re-run later with `RUN_POST_TRIM_MQC=0` to skip the extra MultiQC, and `RUN_TRIM_PASS2=0` if pass-2 is no longer needed.

## `samples.tsv`

Optional batch metadata for the combined project lives next to `data/` in that project directory (`fastq_prefix`, `batch`, `batch_rep_label`, etc.).
