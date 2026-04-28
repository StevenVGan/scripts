# PRO-seq / GRO-seq upstream pipeline

Fork of [`../csRNA`](../csRNA) adapted for PRO-seq chemistry. The csRNA fork already handles stranded small-insert libraries correctly; this fork adds three PRO-seq-specific changes:

1. **Strand flip** (`PROSEQ_FLIP_STRAND=1`) — PRO-seq R1 is the reverse complement of the nascent RNA (Kwak et al. 2013, *Science*). `makeTagDirectory -flip` + swapped `--filterRNAstrand` labels make `_fwd.bw` and all HOMER tag-based tools reflect the RNA strand.
2. **Feature calling** — `findPeaks -style groseq` (via `HOMER_STYLE`) replaces `findcsRNATSS.pl`; purpose-built for GRO/PRO-seq nascent transcripts.
3. **Pausing index + divergent eRNA** (step 4.3) — per-gene PI = promoter density / gene-body density, plus bidirectional pairs from groseq transcripts. Replaces csRNA's `4.3_tss_csrna.sh`.

## Workflow

```
Raw FASTQs (data/)  →  1_trim_qc  →  1.1_polya_trim (optional)  →  2_bowtie2  →  3_homer_tags  →  4.1 / 4.2 / 4.3  →  5_qc
```

**Step 1.1** (opt-in) runs cutadapt on already-trimmed reads to strip poly-A runs at 3′ ends and poly-T runs at 5′ ends of R1 (and R2 for PE). PRO-seq R1 is the reverse complement of the nascent RNA, so mRNA poly-A carryover appears as poly-T at R1's 5′ end; genuine poly-A in the reads comes from A-tracts or sense-oriented contaminants. Enable per-project with `RUN_POLYA_TRIM=1` when FastQC flags either pattern.

## Key variables (`0_config.sh` or env)

| Variable | Default | Meaning |
|---|---|---|
| `SE` | `0` | `1` for single-end projects |
| `RUN_POLYA_TRIM` | `0` | Enable step 1.1 poly-A / poly-T trim after step 1 |
| `POLYA_MIN_STRETCH` | `10` | Minimum A/T run length (cutadapt `-O`) for step 1.1 |
| `POLYA_ERROR_RATE` | `0.1` | Cutadapt `-e` for step 1.1 |
| `PROSEQ_FLIP_STRAND` | `1` | Flip read strand → RNA strand (HOMER `-flip`, bigWig label swap) |
| `HOMER_STYLE` | `groseq` | Passed to `findPeaks -style` (step 4.2) |
| `HOMER_SS_PE` | `1` | `-sspe` for PE tag dirs (forced 0 when `SE=1`) |
| `BT2_EXTRA` | empty | Stranded PE orientation: `--fr` (TruSeq PRO-seq) or `--rf` |
| `TRIM_MIN_LENGTH` | `20` | Short-insert floor; inserts typically 20–65 nt |
| `SMALL_RNA_5P_ADAPTER` | TruSeq SmallRNA | R2 3′ readthrough; blank to disable |
| `STRAND_BIGWIG` | `1` | Emit `_fwd.bw` and `_rev.bw` |
| `COMBINED_BIGWIG` | `1` | Also emit unstranded `${sample}.bw` |
| `BAMCOV_IGNORE_DUP` | `0` | PRO-seq dynamic range resembles RNA-seq — no dedup for coverage |
| `GENE_ANNOT_BED` | `~/work/ref/annot/hg38/hg38_annot_genome.bed` | BED(12) of transcripts (first 6 cols used) for step 4.3 |
| `TSS_WINDOW` | `250` | Promoter half-width (bp around TSS) for PI numerator |
| `GB_START` | `500` | Gene body offset downstream of TSS |
| `GB_MIN_LEN` | `1000` | Minimum gene length to include in PI |
| `DIVERGENT_WINDOW` | `1000` | Max bp upstream for antisense pair detection |
| `RUN_PEAK_MACS3` | `0` | Rarely appropriate for nascent RNA |
| `RUN_PEAK_HOMER` | `1` | `findPeaks -style groseq` |
| `RUN_PAUSING` | `1` | Step 4.3 (pausing index + divergent pairs) |

## Why PRO-seq needs the strand flip

PRO-seq library prep (Kwak et al. 2013) ligates the **5′ Illumina adapter to the RNA's 3′ end** (the biotin-NTP / Pol II position) and the **3′ Illumina adapter to the RNA's 5′ end**. Reverse transcription then primes from the 3′ adapter, producing a cDNA that is the reverse complement of the RNA. R1 reads this cDNA, so R1 aligns to the **opposite** genomic strand from the gene's sense strand. Without flipping, `${sample}_fwd.bw` shows signal for genes on the − strand. The same logic applies to HOMER tag strands, so `findPeaks -style groseq` also expects RNA-strand tags.

## Feature calling

`findPeaks -style groseq` outputs contiguous nascent transcript intervals from stranded tag directories and handles background without a separate input library (PRO-seq has no natural input — nuclear run-on produces its own signal). If `PEAKCALL_GROUPS_FILE` does not exist, step 4.2 auto-creates a no-control groups file with one row per BAM.

Output: `${HOMER_PEAK_DIR}/<name>.annotatePeaks.txt` (blacklist-filtered). Intermediate `*_transcripts.txt` is collapsed to `*_transcripts_stats.txt` (headers only).

## Pausing index (step 4.3)

```
PI = (promoter_count / promoter_length) / (genebody_count / genebody_length)
   promoter = TSS ± TSS_WINDOW
   genebody = TSS+GB_START .. TTS  (only genes ≥ GB_MIN_LEN)
```

Counts are from `bedtools multicov` with strand-aware flag (`-S` when `PROSEQ_FLIP_STRAND=1`, else `-s`). Outputs per sample:

- `${sample}_pausing_index.tsv` — per-gene table (gene_id, strand, lengths, counts, densities, PI)
- `${sample}_divergent_pairs.tsv` — closest antisense pair per + strand transcript (from 4.2 output) within `DIVERGENT_WINDOW`
- `pausing_index_summary_mqc.tsv` — per-sample N, median, Q1, Q3 for MultiQC custom content

## Running

From a project's `script/` folder:

```bash
./run_all.sh                                     # full pipeline
SE=1 ./run_all.sh                                # single-end
RUN_TRIM=0 RUN_BOWTIE2=0 RUN_PAUSING=1 ./run_all.sh   # redo pausing only
```

## Strand sanity check (IGV)

- `GAPDH` (chr12, + strand) → signal in `${sample}_fwd.bw`
- `ACTB` (chr7, − strand) → signal in `${sample}_rev.bw`

If these are reversed, set `PROSEQ_FLIP_STRAND=0` and re-run steps 2 + 3 (everything downstream will pick up the fresh tag dirs and bigWigs).

## Out of scope

- dREG enhancer calling (GPU/Python deps)
- Differential transcription analysis (DESeq2 on gene-body counts) — belongs in per-project `analysis/` folder
- grohmm / NRSA R packages
