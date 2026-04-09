# Prep tools (upstream / project setup)

Scripts for downloading raw FASTQs, merging Illumina lanes, and linking into a CUTRUN (or other) project **`data/`** tree before **`1_trim_qc.sh`**.

| Script | Purpose |
|--------|---------|
| **download_fastq.sh** | `wget` from IGM FTP + `md5sum -c`; edit `FTP_BASE`, `DEST_DIR`, `FILE_PREFIX` |
| **link_fastq.sh** | Symlink from `RAW_DIR` using `prefix`→`newname` map (`*_R1_001` / `*_R2_001` or SE `.fastq.gz`) |
| **merge_lanes_inplace.sh** | Concatenate lanes in **`IGM_DIR`**, write `${sample}_R{1,2}.fastq.gz`, **delete** lane files; needs `IGM_DIR` + `MAP_FILE` |
| **link_merged_fastqs.sh** | Symlink merged `${sample}_R{1,2}.fastq.gz` into **`DEST_DIR`**; set **`MERGED_FASTQ_DIR`** or **`IGM_FASTQ_DIR`** (same folder); needs **`MAP_FILE`** |
| **download_geo_fastq_ena.sh** | ENA HTTPS `curl` for SRR list (SE-only per SRR); parallel **`DOWNLOAD_JOBS`** |

Typical **multi-lane IGM** flow:

```bash
# 1) Edit config in download_fastq.sh, then:
./download_fastq.sh

# 2) Merge + delete lanes (destructive in IGM_DIR)
IGM_DIR=~/work/raw_seq/RUN_ID MAP_FILE=~/work/seq/CUTRUN/proj/igm_to_sample.tsv ./merge_lanes_inplace.sh

# 3) Link into project data/
MERGED_FASTQ_DIR=~/work/raw_seq/RUN_ID DEST_DIR=~/work/seq/CUTRUN/proj/data \
  MAP_FILE=~/work/seq/CUTRUN/proj/igm_to_sample.tsv ./link_merged_fastqs.sh
```

**Single-lane** (or already one file per mate): skip merge; use **`link_fastq.sh`** with `link_sample.tsv`.

**GEO / ENA (SRR, single-end):** edit config in **`download_geo_fastq_ena.sh`**, then run from this directory (or call with absolute path: `bash /path/to/prep/download_geo_fastq_ena.sh`).

**Paths:** All scripts in this folder are path-agnostic; use absolute paths in **`MAP_FILE`** / env vars if you do not `cd` here.

Parent index: [../README.md](../README.md).
