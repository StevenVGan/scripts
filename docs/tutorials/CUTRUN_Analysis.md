---
layout: default
title: CUT&RUN Analysis Tutorial
---

# CUT&RUN Analysis Tutorial: Mapping Histone Modifications and Transcription Factor Binding Sites

**Author:** Steven Gan

**Date:** February 2026

**Objective:** Learn how to analyze CUT&RUN sequencing data using command-line tools on an AWS EC2 instance. This tutorial covers data acquisition, quality control, read alignment, peak calling, and visualization.

---

## Background

Cleavage Under Targets and Release Using Nuclease (CUT&RUN) is a powerful technique used to identify genome-wide binding sites of transcription factors, histone modifications, and other DNA-associated proteins. Unlike traditional ChIP-seq, CUT&RUN uses an antibody-targeted nuclease to cleave DNA adjacent to the protein of interest, resulting in lower background and requiring fewer cells.

In this tutorial, you will analyze CUT&RUN data from the paper ["Different modes of engagement with the nucleosome acidic patch yield distinct functional outcomes"](https://www.biorxiv.org/content/10.64898/2026.01.14.699602v1). The researchers performed CUT&RUN for Dot5, a Dot1-like protein in *Saccharomyces cerevisiae* (budding yeast). They also included an H3K4me3 positive control, which is a histone modification associated with active promoters.

The study found that "CUT&RUN failed to identify any regions of Dot5-FLAG enrichment" despite showing that Dot5 is associated with chromatin via cell fractionation. This makes it an interesting dataset for learning the analysis workflow.

**The CUT&RUN workflow involves several key steps:**

- Download sequencing data from NCBI SRA
- Read trimming to remove adapter sequences and low-quality bases  
- Quality control of trimmed sequencing reads
- Alignment of paired-end reads to a reference genome
- Peak calling to identify regions of enrichment
- Quality assessment with MultiQC
- Peak annotation with HOMER
- Visualization of binding patterns with heatmaps

**Important note:** We will be working with **paired-end sequencing data**, which means each DNA fragment is sequenced from both ends, giving us two reads (R1 and R2) per fragment.

---

## AWS EC2 Instance Setup

**Recommended specifications:**

- **AMI:** Ubuntu Server 24.04 LTS (HVM), SSD Volume Type
- **Instance Type:** m5.xlarge (4 vCPUs, 16 GiB memory)
- **Storage:** 60 GB General Purpose SSD (gp3)
- **Security Group:** Your class security group
- **Key Pair:** Use your existing key pair
- **Tag:** Add `terminable=true` for cleanup

**Important:** For Ubuntu instances, the username is **`ubuntu`**, not `ec2-user`!

---

## STEP 1: Setting Up the Conda Environment

**Preliminary Note:** Before everything starts, it is a good idea to keep all the commands that you run in a text file for future reference. You can create a text file using the command `touch commands.txt` and then edit it easily.

### Install Miniconda

Download and install Miniconda:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Answer "yes" and press enter to all prompts during installation.

Reload your shell configuration or just close and reopen your terminal:
```bash
source ~/.bashrc
```

---

### Configure Conda Channels

Add channels in the correct order:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

**Command breakdown:**

- `conda config --add channels <name>`: Adds a package repository to conda
- Order matters: conda-forge > bioconda > defaults (priority from top to bottom)

Check that channels are in the correct order:
```bash
conda config --show channels
```

Expected output:
```
channels:
  - conda-forge
  - bioconda
  - defaults
```

---

### Accept Terms of Service

Run these commands to accept the required terms (if warning messages appear during package installation):
```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

---

### Create Environment with Required Packages

```bash
conda create -n cutrun fastqc trim-galore bowtie2 macs3 samtools bedtools deeptools sra-tools multiqc homer -y
```

**Command breakdown:**

- `conda create -n cutrun`: Creates a new isolated environment named "cutrun"
- `fastqc`: Quality control tool for sequencing data
- `trim-galore`: Wrapper for Cutadapt and FastQC for adapter trimming
- `bowtie2`: Fast read aligner using Burrows-Wheeler transform
- `macs3`: Model-based Analysis of ChIP-Seq version 3 (peak caller)
- `samtools`: Suite of tools for manipulating SAM/BAM alignment files
- `bedtools`: Toolkit for genome arithmetic and interval operations
- `deeptools`: Tools for visualizing and analyzing deep-sequencing data
- `sra-tools`: NCBI tools for downloading data from Sequence Read Archive
- `multiqc`: Aggregates results from bioinformatics analyses
- `homer`: Hypergeometric Optimization of Motif EnRichment (annotation tool)
- `-y`: Automatically answers "yes" to installation prompts

**Note:** Installation takes 5-10 minutes.

---

### Activate the Environment

```bash
conda activate cutrun
```

Your prompt should now show `(cutrun)` at the beginning.

Check installed package versions:
```bash
fastqc --version
bowtie2 --version
```

---

### Create Directory Structure


Now we will create the directory structure for the analysis. We will have a main project directory called `~/cutrun/` with subdirectories for raw data, trimmed data, alignment results, coverage tracks, peak calling results, and matrices for heatmaps. We will also create a reference genome directory `~/ref/`. The reference directory is shared across projects, while the `cutrun` directory contains all files specific to this analysis.

```bash
mkdir ~/ref
mkdir ~/cutrun
cd ~/cutrun
mkdir data cleandata align peaks matrix tracks
```

**Directory structure explanation:**
```
~/ref/              # Reference genomes (shared across projects)
~/cutrun/           # Main project directory
  ├── data/         # Raw FASTQ files from SRA
  ├── cleandata/    # Trimmed/filtered FASTQ files
  ├── align/        # BAM alignment files
  ├── peaks/        # Peak calling results
  ├── matrix/       # Matrices for heatmaps
  └── tracks/       # BigWig coverage tracks
```

**Command breakdown:**

- `mkdir ~/ref`: Creates a reference genome directory in your home directory
- `mkdir ~/cutrun`: Creates the main project directory
- `cd ~/cutrun`: Changes to the project directory
- `mkdir data cleandata align peaks matrix tracks`: Creates all subdirectories at once

---

## STEP 2: Data Acquisition

### About the Dataset

Data are available from GEO under accession number **GSE313255**

There are 4 samples:

- **GSM9364728** - IgG Control - (SRR36395843)
- **GSM9364729** - anti-H3K4me3 - (SRR36395842)
- **GSM9364730** - anti-FLAG Control - (SRR36395841)
- **GSM9364731** - DOT5-FLAG - (SRR36395840)

All samples are paired-end sequencing data.

---

### Download Data from SRA

Navigate to the data directory:
```bash
cd ~/cutrun/data/
```

Download all four samples using `fasterq-dump`:
```bash
fasterq-dump SRR36395843 SRR36395842 SRR36395841 SRR36395840 --split-files -O .
```

**Command breakdown:**

- `fasterq-dump`: Fast FASTQ download tool from SRA toolkit
- `SRR36395843 SRR36395842 SRR36395841 SRR36395840`: Four SRA accession numbers (space-separated)
- `--split-files`: Separates paired-end reads into _1.fastq and _2.fastq files
- `-O .`: Outputs files to the current directory


Check the downloaded files:
```bash
ls -lh
```

---

### Rename Files Using Arrays and For Loop

The SRR accession numbers are not descriptive. Let's rename them to meaningful sample names. 

<!-- SRR36395843 -> IgG_Control
SRR36395842 -> H3K4me3
SRR36395841 -> FLAG_Control
SRR36395840 -> DOT5-FLAG -->

**Renaming Scheme:**

- SRR36395843 -> IgG_Control
- SRR36395842 -> H3K4me3
- SRR36395841 -> FLAG_Control
- SRR36395840 -> DOT5-FLAG

We can do it one by one:
```bash
mv SRR36395843_1.fastq IgG_Control_R1.fastq
mv SRR36395843_2.fastq IgG_Control_R2.fastq
mv SRR36395842_1.fastq H3K4me3_R1.fastq
# ...
```

But it's more efficient to use arrays and a for loop:

**Step 1: Define the files and names in arrays**
```bash
FILES=("SRR36395843" "SRR36395842" "SRR36395841" "SRR36395840")
NAMES=("IgG_Control" "H3K4me3" "FLAG_Control" "DOT5-FLAG")
```

**Command breakdown:**

- `FILES=(...)`: Creates a variable called FILES that contains all the SRR accession numbers
- `NAMES=(...)`: Creates a variable called NAMES that contains the corresponding sample names
- Note: The `=` sign has **no spaces** around it!
- Customarily, variables are in UPPERCASE letters to distinguish them from commands
- To print a variable: `echo ${FILES[0]}` will print `SRR36395843`
- Note how `$` and `{}` are used to access the variable, and `[]` accesses the array index (starting from 0).

**Step 2: Loop through the arrays and rename files**

For loops have the syntax:

```bash
for VARIABLE in LIST; do <commands>; done
```

where `VARIABLE` takes on each value in `LIST` one at a time. `;` is used to separate commands on the same line, or you can put them on separate lines.

Now loop through the indices of the arrays to rename files:
```bash
for i in ${!FILES[@]}; do
    mv ${FILES[$i]}_1.fastq ${NAMES[$i]}_R1.fastq
    mv ${FILES[$i]}_2.fastq ${NAMES[$i]}_R2.fastq
done
```

**Command breakdown:**

- `${!FILES[@]}`: The `@` expands all the elements in the FILES array, and the `!` extracts their indices (0, 1, 2, 3).
- `i`: Loop variable that takes on each value from `${!FILES[@]}`, which are the indices 0, 1, 2, 3
- `${FILES[$i]}`: Accesses the element at index i (e.g., `${FILES[0]}` = "SRR36395843")
- `mv ${FILES[$i]}_1.fastq ${NAMES[$i]}_R1.fastq`: Renames the R1 file

**Step 3: Compress the raw FASTQ files to save space**
```bash
gzip *.fastq
```

**Command breakdown:**

- `gzip`: Command to compress files (creates .gz files)
- `*`: Wildcard that matches all files ending with .fastq
- This takes ~ 7 minutes for all files
- Most bioinformatics tools can read gzipped files directly, so no need to unzip for analysis
- Check the file sizes: `ls -lh` - files should be much smaller now

Verify the renaming and compression:
```bash
ls -lh *.fastq.gz
```


---

## STEP 3: Trimming and Post-trimming Quality Control

### Understanding Trim Galore

Trim Galore is a wrapper around Cutadapt and FastQC that:

- Automatically detects and removes adapter sequences
- Trims low-quality bases from the 3' end
- Removes reads that are too short after trimming
- Runs FastQC on the trimmed files automatically

**Note:** We perform trimming BEFORE standalone quality control because trim_galore runs FastQC automatically, and we want to assess the quality of the data we'll actually use for alignment.

---

### Trim All Samples Using a For Loop

Navigate to the work directory:
```bash
cd ~/cutrun/
```

Trim_galore has the following syntax for paired-end data:
```bash
trim_galore --paired <R1.fastq.gz> <R2.fastq.gz> -o <output_directory>
```

It also accepts the arguments such as `--stringency <value>, --length <value>, --quality <value>` to customize trimming parameters. For this analysis, we will use stringency 5, length 30, quality 28, which is more stringent than default.


Use a for loop to trim all samples:
```bash
SAMPLES=("IgG_Control" "H3K4me3" "FLAG_Control" "DOT5-FLAG")

for sample in ${SAMPLES[@]}; do
  trim_galore --paired \
    --stringency 5 --length 30 --quality 28 \
    ./data/${sample}_R1.fastq.gz ./data/${sample}_R2.fastq.gz \
    -o ./cleandata/
done
```

**Command breakdown:**

- `SAMPLES=(...)`: Creates an array of sample names
- `for sample in ${SAMPLES[@]}; do ... done`: Loops through each sample name
- `\` is used to split one long commands into multiple lines for readability. This is different from `;`, which separates multiple commands to be run sequentially. It is important to note that the space before `\` is required, as it represents the space which exists if the command were written on a single line.
- `trim_galore --paired`: Runs Trim Galore in paired-end mode
- `-o ./cleandata/`: Specifies the output directory for trimmed files.


---


### Understanding the Output Files

After trimming, you'll have these files for each sample:
- `*_R1_val_1.fq`: Trimmed R1 reads
- `*_R2_val_2.fq`: Trimmed R2 reads
- `*_R1.fastq.gz_trimming_report.txt`: Trimming statistics for R1
- `*_R2.fastq.gz_trimming_report.txt`: Trimming statistics for R2

**Why "_val_1" and "_val_2"?** Trim Galore adds these suffixes to indicate "validated" reads that passed quality filters.

Most reads were retained after trimming, indicating that the data quality was generally good!


---

### Post-trimming Quality Control with FastQC

Now that we have trimmed the reads, we can perform quality control using FastQC to assess the quality of the trimmed reads. FastQC provides a quick overview of the quality of sequencing data, including per-base quality scores, GC content, sequence duplication levels, and overrepresented sequences.

We will first inspect the quality of the trimmed reads using FastQC. Create a subdirectory for FastQC output in `./cutrun/cleandata/`:
```bash
mkdir ./cleandata/fastqc/
```

Run FastQC on all gzipped FASTQ files:
```bash
fastqc ./cleandata/*.gz -o ./cleandata/fastqc/
```

**Command breakdown:**

- `./cleandata/*.gz`: Input all gzipped FASTQ files in the cleandata directory
- `-o ./cleandata/fastqc/`: Output directory for FastQC results

We can either inspect the FastQC results locally by downloading the HTML files, or use an extension in VSCode to view HTML files directly on the EC2 instance. Install the "HTML Preview" extension in VSCode, then right-click on the HTML file and select "Preview".


---

## STEP 4: Reference Genome Preparation

We will align the reads to the *S. cerevisiae* reference genome (SacCer3). Usually in a well managed lab environment, the reference genomes are stored in a shared directory, e.g. `~/ref/`. For the practice here, we will download the reference genome and build the Bowtie2 index ourselves.

### Download Yeast Reference Genome

Creeat a subdirectory for Genome index and navigate to the reference directory:
```bash
mkdir ~/ref/sacCer3
cd ~/ref/sacCer3
```

Download the *S. cerevisiae* reference genome:
```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
```

**Command breakdown:**

- `wget`: Downloads files from the web
- URL: Direct link to the yeast reference genome from Saccharomyces Genome Database (SGD)
- `gunzip`: Decompresses .gz files, opposite of `gzip`.
- Output file will have the same name without .gz extension


### Build Bowtie2 Index

Bowtie2 requires an index of the reference genome for efficient alignment:

```bash
bowtie2-build sacCer3.fa sacCer3
```

**Command breakdown:**

- `bowtie2-build`: Builds a Bowtie2 index using the Burrows-Wheeler transform
- `sacCer3.fa`: Input reference genome in FASTA format
- `sacCer3`: Prefix/basename for the index files

This creates 6 index files with extensions `.bt2`:

- `sacCer3.1.bt2`, `sacCer3.2.bt2` (forward index)
- `sacCer3.3.bt2`, `sacCer3.4.bt2` (reverse index)
- `sacCer3.rev.1.bt2`, `sacCer3.rev.2.bt2` (reverse complement index)

**Note:** For the small yeast genome (~12 Mb), this takes about 1 minute. For human genome (~3 Gb), it would take ~ 2 hours! This is why reference genomes are usually pre-indexed and shared in a lab environment.


---

## STEP 5: Read Alignment with Bowtie2

### Understanding Alignment

Alignment maps each sequencing read to its most likely genomic location. For paired-end data:

- R1 and R2 reads from the same fragment should map to opposite strands
- They should be within a reasonable distance (insert size)
- "Concordant" pairs have both mates mapping correctly

---

### Align All Samples Using a For Loop

Navigate to the alignment directory:
```bash
cd ~/cutrun/
```

Bowtie2 command for paired-end alignment:
```bash
bowtie2 -x <index_basename> \
    -1 <R1.fastq> -2 <R2.fastq> \ 
    -S <output.sam> \
    --very-sensitive \
    -p <num_threads>
```

**Command breakdown:**

- `-x`: Basename of the Bowtie2 index files
- `-1`: Input R1 FASTQ file
- `-2`: Input R2 FASTQ file
- `-S`: Output SAM file
- `--very-sensitive`: Preset for high sensitivity alignment (slower but better)
- `-p`: Number of CPU threads to use (default is 1). This server has 4 vCPUs, so we can use 4 threads. But to avoid overloading the server, we will use 3 threads instead.


Align all samples with a for loop (we have defined the SAMPLES array earlier, if you are unsure, you can print the SAMPLES array using `echo ${SAMPLES[@]}`):
```bash

GENOME=~/ref/sacCer3/sacCer3

for sample in ${SAMPLES[@]}; do
  bowtie2 -x $GENOME \
    -1 ./cleandata/${sample}_R1_val_1.fq.gz \
    -2 ./cleandata/${sample}_R2_val_2.fq.gz \
    -S ./align/${sample}.sam \
    --very-sensitive \
    -p 3
done
```

**Note:** I put genome path in a variable `$GENOME` for easier reuse, since this is usually a long path...

Now it is always scary at this step because alignment can take a long time and `bowtie2` does not create any output files until it is finished. To make sure it is running, we can open another terminal and use `top` or `htop` command to monitor CPU usage (type `q` to exit `top` or `htop`).

You can also check the size of the SAM files in ./align/ using `ls -lh ./align/` to see if they are being created and growing in size.


---

## STEP 6: SAM/BAM Processing

### Understanding SAM and BAM Formats

- **SAM** (Sequence Alignment/Map): Human-readable text format containing alignment information
- **BAM**: Compressed binary version of SAM (saves ~80% disk space)
- **Sorting**: Organizes reads by genomic coordinate (required for most downstream tools)
- **Indexing**: Creates a .bai index file for rapid random access to specific regions

---

### Convert, Sort, and Index BAM Files

Process all samples with a for loop:
```bash
for SAMPLE in ${SAMPLES[@]}; do
    samtools flagstat ./align/${SAMPLE}.sam > ./align/${SAMPLE}_flagstat.txt
    samtools view -bS ./align/${SAMPLE}.sam -o ./align/${SAMPLE}.bam
    samtools sort ./align/${SAMPLE}.bam -o ./align/${SAMPLE}_sorted.bam
    samtools index ./align/${SAMPLE}_sorted.bam
done
```

**Command breakdown:**

- `samtools flagstat <input.sam> > <output_flagstat.txt>`: Generates alignment statistics
  - `flagstat`: Command to compute alignment statistics
  - `>`: Redirects output to a text file, as the default is to print to terminal. Other programs usually use `-o` for output files.

- `samtools view -bS <input.sam> -o <output.bam>`: Converts SAM to BAM
  - `view`: Command to view/convert alignment files
  - `-bS`: Input is SAM, output is BAM

- `samtools sort <input.bam> -o <output_sorted.bam>`: Sorts BAM file by genomic coordinates
  - `sort`: Command to sort BAM files

- `samtools index <input_sorted.bam>`: Creates index file for sorted BAM
  - `index`: Command to index BAM files
  - Generates a .bai file alongside the sorted BAM. Indexing allows rapid access to specific genomic regions.


---

### Clean Up Intermediate Files

Check whether all sorted BAM files are created:
```bash
ls -lh ./align/
```

Remove SAM files and unsorted BAM files to save disk space:

```bash
for SAMPLE in ${SAMPLES[@]}; do
    rm ./align/${SAMPLE}.sam
    rm ./align/${SAMPLE}.bam
done
```

**Note:** DO NOT USE `rm ./align/*.sam` or `rm ./align/*.bam` as it will also delete the sorted BAM files as they also end with `.bam`!

---

### Examine Alignment Statistics

Check the flagstat output files:
```bash
for SAMPLE in ${SAMPLES[@]}; do
    echo "=== ${SAMPLE} ==="
    cat ./align/${SAMPLE}_flagstat.txt
    echo ""
done
```

It seems that H3K4me3 has the highest mapping rate (92.18%), while the other samples have lower mapping rates (~65%).

---

## STEP 7: Creating BigWig Coverage Tracks

### Understanding BigWig Files

BigWig files store genome-wide coverage data in a compressed, indexed format:
- Much smaller than BAM files
- Can be directly loaded into genome browsers (IGV, UCSC Genome Browser)
- Represent normalized coverage values across the genome
- Useful for visualization and comparison between samples

---

### Generate BigWig Files for All Samples

DeepTools command syntax to generate bigWig files:
```bash
bamCoverage \
    -b <input.bam> \
    -o <output.bw> \
    --binSize <value> \
    --normalizeUsing <method> \
    --ignoreDuplicates \
    -p <num_threads>
```

**Command breakdown:**

- `-b`: Input BAM file
- `-o`: Output BigWig file
- `--binSize`: Size of bins to calculate coverage (smaller = higher resolution, larger = smaller file)
- `--normalizeUsing`: Normalization method (e.g., RPKM, CPM, BPM, RPGC)
- `--ignoreDuplicates`: Ignores PCR duplicates in coverage calculation
- `-p`: Number of CPU threads to use

Create BigWig files using bamCoverage from deepTools:
```bash
for sample in ${SAMPLES[@]};
do
  bamCoverage \
    -b ./align/${sample}_sorted.bam \
    -o ./tracks/${sample}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --ignoreDuplicates \
    -p 3
done
```

Download the bigWig files to your local computer for visualization in IGV or UCSC Genome Browser ...

From the tracks, we can tell that the H3K4me3 works well and shows strong enrichment at promoter regions, while the DOT5-FLAG sample does not show any clear enrichment, consistent with the original paper's claim.



---

## STEP 8: Peak Calling with MACS3

### Understanding Peak Calling

Peak calling identifies genomic regions with significant enrichment in the treatment sample compared to the control:
- Compares **treatment** (ChIP/CUT&RUN sample) to **control** (IgG/Input)
- Uses statistical models to determine significance
- Outputs enriched regions (peaks) in various formats

For our samples:
- **H3K4me3** (treatment) vs. **IgG_Control** (control)
- **DOT5-FLAG** (treatment) vs. **FLAG_Control** (control)

---

### Call Peaks for H3K4me3

MACS3 command syntax for paired-end data:
```bash
macs3 callpeak \ 
  -t <treatment.bam> \
  -c <control.bam> \
  -f BAMPE \
  -g <genome_size> \
  -n <output_prefix> \
  --outdir <output_directory> \
  -q <qvalue_threshold>
```

**Command Breakdown:**

- `-t`: Treatment BAM file
- `-c`: Control BAM file
- `-f`: Input format (BAMPE for paired-end BAM files)
- `-g`: Effective genome size (1.2e7 for yeast)
- `-n`: Output prefix for peak files
- `--outdir`: Output directory
- `-q`: Q-value threshold for peak calling (0.01)


For H3K4me3 sample, we will use IgG_Control as the control sample.

```bash
macs3 callpeak \
  -t ./align/H3K4me3_sorted.bam \
  -c ./align/IgG_Control_sorted.bam \
  -f BAMPE \
  -g 1.2e7 \
  -n H3K4me3 \
  --outdir ./peaks/ \
  -q 0.01
```

For DOT5-FLAG sample, we will use FLAG_Control as the control sample.

```bash
macs3 callpeak \
  -t ./align/DOT5-FLAG_sorted.bam \
  -c ./align/FLAG_Control_sorted.bam \
  -f BAMPE \
  -g 1.2e7 \
  -n DOT5-FLAG \
  --outdir ./peaks/ \
  -q 0.01
```

### Examine Peak Calling Output

We can list the number of peaks called for .narrowPeak files:
```bash
wc -l ./peaks/*.narrowPeak
```

Expected output:
```
     1 ./peaks/DOT5-FLAG_peaks.narrowPeak
  3031 ./peaks/H3K4me3_peaks.narrowPeak
```

**Interpretation:**
- H3K4me3 shows 3031 peaks (strong enrichment, as expected for this histone mark)
- DOT5-FLAG shows only 1 peak (essentially no enrichment, consistent with the paper's findings)

---

## STEP 9: Quality Assessment with MultiQC

### Understanding MultiQC

MultiQC aggregates results from multiple bioinformatics tools into a single interactive HTML report:

- Combines FastQC, Samtools, Bowtie2 outputs
- Allows easy comparison across samples
- Helps identify outliers or quality issues
- Generates summary statistics tables and plots

---

### Configure MultiQC

Create a configuration file `multiqc_config.yml` to properly merge related samples:

Add the following content:
```yaml
table_sample_merge:
  "Raw R1": "_R1"
  "Raw R2": "_R2"
  "Trimmed R1": "_R1_val_1"
  "Trimmed R2": "_R2_val_2"
  "Alignment": "_flagstat"
```

**Configuration explanation:**
- Groups FastQC reports by read type instead of treating each file separately.
- Merges flagstat reports under "Alignment"
- Makes the MultiQC report cleaner and easier to interpret

Save and exit.

---

### Run MultiQC

```bash
multiqc . -o ./multiqc/ -c multiqc_config.yml
```

**Command breakdown:**

- `.`: Input directory (current directory, multiqc will search all subdirectories)
- `-o`: Output directory
- `-c`: Specifies the custom configuration file

The MultiQC report will be saved in `./multiqc/multiqc_report.html`. Download this HTML file to your local computer to view the interactive report.

---

## STEP 10: Peak Annotation with HOMER

### Understanding Peak Annotation

Peak annotation assigns genomic features to each peak:
- **Promoter-TSS**: Within ~1 kb of transcription start sites
- **Exon**: Within protein-coding exons
- **Intron**: Within introns
- **TTS**: Near transcription termination sites
- **Intergenic**: Between genes

This helps interpret the biological function of the binding sites.

---

We can use HOMER to annotate the called peaks to genomic features such as promoters, exons, introns, intergenic regions, etc. This helps us understand the distribution of the peaks across different genomic regions.

In order to use HOMER for peak annotation, we need a gene annotation file in GTF format for *S. cerevisiae*. We can download it from Ensembl or UCSC.


### Download Gene Annotation



Download the GTF annotation file for yeast:
```bash
cd ~/ref/
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
gunzip sacCer3.ensGene.gtf.gz
```

**Command breakdown:**
- Downloads gene annotations in GTF format from UCSC Genome Browser
- `gunzip`: Decompresses the file

---

### Configure HOMER

Configure HOMER to recognize the sacCer3 genome assembly:
```bash
perl ~/miniconda3/envs/cutrun/share/homer/configureHomer.pl -install sacCer3
```

**Command breakdown:**

- If HOMER `annotatePeaks.pl` cannot find the genome, it will prompt you to run this command to install necessary data files for the specified genome assembly (sacCer3 in this case).
- This only needs to be done once per genome assembly.

---

### Annotate H3K4me3 Peaks

annotatePeaks.pl command syntax:
```bash
annotatePeaks.pl <peak_file> <genome> -gtf <annotation_file> > <output_file>
```

**Command breakdown:**
- `<peak_file>`: Input peak file in BED or narrowPeak format
- `<genome>`: Genome assembly name (e.g., sacCer3)
- `-gtf <annotation_file>`: GTF file for gene annotations
- `>`: Redirects output to a text file
- `<output_file>`: Output file containing annotated peaks

Annotate the H3K4me3 peaks:
```bash
cd ~/cutrun/
annotatePeaks.pl ./peaks/H3K4me3_peaks.narrowPeak sacCer3 -gtf ~/ref/sacCer3.ensGene.gtf > ./peaks/H3K4me3_peaks_annotated.txt
```

The annotated file contains detailed information for each peak:
- Peak coordinates
- Nearest gene
- Distance to TSS
- Genomic feature annotation (promoter, exon, etc.)

For practice, let's print the Percentage of each genomic feature, which is in column 8. We want to exclude the "(XXX)" part from the annotation.

Extract and count genomic feature annotations:
```bash
cut -f8 ./peaks/H3K4me3_peaks_annotated.txt | sed 's/ (.*)//' | sort | uniq -c
```

**Command breakdown:**

- `cut -f8`: Extracts column 8 (Annotation column)
- `sed 's/ (.*)//'`: Removes everything after " (" using substitution
  - `s/pattern/replacement/`: Substitute pattern with replacement
  - ` (.*)`: Space followed by opening parenthesis and any characters until closing parenthesis
  - This removes distance information like "(123 bp)"
- `sort`: Sorts annotations alphabetically
- `uniq -c`: Counts unique entries and displays count

Expected output:
```
      1 Annotation
      6 Intergenic
    620 TTS
    134 exon
      3 intron
   2268 promoter-TSS
```

**Interpretation:**

- Annotation is the header line, it is counted one time.
- Most H3K4me3 peaks (2268/3031 = 75%) are at promoter-TSS regions
- This is consistent with H3K4me3's known role as a mark of active promoters
- Some peaks at TTS (transcription termination sites) and within gene bodies (exons)

---

## STEP 11: Visualization with Heatmaps

### Understanding Heatmaps

Heatmaps display signal intensity across many genomic regions:

- **Rows**: Individual genomic regions (e.g., peaks or genes)
- **Columns**: Genomic positions (e.g., ±2 kb from peak center)
- **Color intensity**: Signal strength (read coverage)

This visualization reveals patterns across thousands of regions simultaneously.

---

### Compute Matrix for Peak Regions

First, we will create a matrix of signal values around the peaks using `computeMatrix` tool from deepTools.

DeepTools command syntax to compute matrix:
```bash
computeMatrix reference-point \
  -S <bigwig1> <bigwig2> ... \
  -R <regions.bed> \
  -a <upstream_distance> \
  -b <downstream_distance> \
  -o <output_matrix.gz> \
  --referencePoint <point_type> \
  --skipZeros \
  -p <num_threads>
```

**Command breakdown:**

- `-S`: Input BigWig files (signal tracks) - space-separated list
- `-R`: Regions file (BED format) - using peak regions
- `-a`: Distance AFTER the reference point to include
- `-b`: Distance BEFORE the reference point to include
- `-o`: Output matrix file (gzip compressed)
- `--referencePoint`: Type of reference point (e.g., center, TSS, TES)
- `--skipZeros`: Skip regions with no signal in any sample
- `-p`: Number of CPU threads to use


We will plot all four samples together for comparison, and use H3K4me3 peaks as the reference regions. Before that, let's create variables to hold the input BigWig files and labels for easier reuse:

```bash
INPUTSAMPLES="./tracks/IgG_Control.bw ./tracks/H3K4me3.bw ./tracks/FLAG_Control.bw ./tracks/DOT5-FLAG.bw"
LABELS="IgG_Control H3K4me3 FLAG_Control DOT5-FLAG"
```

Now compute the matrix:
```bash
computeMatrix reference-point \
  -S $INPUTSAMPLES \
  -R ./peaks/H3K4me3_peaks.narrowPeak \
  -a 2000 \
  -b 2000 \
  -o ./matrix/cutrun_matrix.gz \
  --referencePoint center \
  --skipZeros \
  -p 3
```

**Note:** We use `$INPUTSAMPLES` variable to hold the list of BigWig files for easier reuse. `$LABELS` variable holds the sample labels, which will be used in plotting. Matrix are computed for ±2 kb around peak centers. 

---

### Plot Heatmap for Peaks


Plot heatmap using `plotHeatmap` from deepTools:
```bash
plotHeatmap \
  -m <input_matrix.gz> \
  -o <output_image.png> \
  --colorMap <color_scheme> \
  --samplesLabel <label1> <label2> ...
```

**Command breakdown:**

- `-m`: Input matrix file (from computeMatrix)
- `-o`: Output image file (PNG format)
- `--colorMap`: Color scheme for heatmap (e.g., Reds, Blues, Viridis)
- `--samplesLabel`: Labels for each sample column (space-separated)

Plot the heatmap for the computed matrix:

```bash
plotHeatmap \
  -m ./matrix/cutrun_matrix.gz \
  -o ./matrix/cutrun_heatmap.png \
  --colorMap Reds \
  --samplesLabel $LABELS
```

Now we can view the heatmap image file ./matrix/cutrun_heatmap.png to see the enrichment patterns of the samples around the H3K4me3 peaks. To futher gain insights to the patterns of H3K4me3 on genes. We can also plot it on TSS of genes using appropriate reference regions. We will download a bed file of TSS regions from UCSC Table Browser for SacCer3, and use it as reference regions to plot heatmap and profile for H3K4me3 sample only. It is available from EPDnew database as well.


---

### Download TSS Annotation and Create TSS Heatmap

Download TSS positions from EPDnew database:
```bash
cd ~/ref/
wget -O sacCer3.EPDnew.TSS.bed https://epd.expasy.org/ftp/epdnew/S_cerevisiae/current/Sc_EPDnew.bed
```

Convert space delimiters to tabs (required for deepTools):
```bash
sed -i 's/ /\t/g' sacCer3.EPDnew.TSS.bed
```

**Command breakdown:**

- `sed -i`: Edit file in-place (modifies the file directly)
- `'s/ /\t/g'`: Substitute (s) spaces with tabs (\t), globally (g = all occurrences)

Navigate back to project directory and compute matrix for TSS regions:
```bash
cd ~/cutrun/
computeMatrix reference-point \
  -S $INPUTSAMPLES \
  -R ~/ref/sacCer3.EPDnew.TSS.bed \
  -a 1000 \
  -b 1000 \
  -o ./matrix/cutrun_TSS.gz \
  --referencePoint TSS \
  --skipZeros \
  -p 3
```

Plot TSS-centered heatmap:
```bash
plotHeatmap \
  -m ./matrix/cutrun_TSS.gz \
  -o ./matrix/cutrun_TSS_heatmap.png \
  --colorMap Reds \
  --samplesLabel $LABELS
```

Now we can see clearly how H3K4me3 is enriched at TSS regions with a very specific pattern that shows the nucleosome-depleted region (NDR) at TSS and well-positioned nucleosomes downstream. This is not observed if plotted at the center of H3K4me3 peaks.

---

## Summary

### Complete Analysis Workflow

In summary, we have completed the CUT&RUN data analysis workflow, including quality control, trimming, alignment, peak calling, annotation, and visualization. This workflow can be adapted for other CUT&RUN datasets with appropriate modifications based on the specific experimental design and biological questions. This tutorial provides a comprehensive guide to analyzing CUT&RUN data using commonly used bioinformatics tools and best practices.

---

### Key Commands Summary

| Tool | Purpose | Key Options |
|------|---------|-------------|
| `fasterq-dump` | Download FASTQ from SRA | `--split-files` (paired-end) |
| `trim_galore` | Trim adapters and low-quality bases | `--paired` `--stringency`, `--length`, `--quality`, `-o` |
| `fastqc` | Quality control of FASTQ files | `-o` (output directory) |
| `bowtie2-build` | Build Bowtie2 index from reference genome | N/A |
| `bowtie2` | Align reads to reference genome | `-x`, `-1`, `-2`, `-S`, `--very-sensitive`, `-p` |
| `samtools view` | Convert SAM to BAM | `-bS`, `-o` |
| `samtools sort` | Sort BAM file by genomic coordinates | `-o` |
| `samtools index` | Index sorted BAM file | N/A |
| `samtools flagstat` | Generate alignment statistics | N/A |
| `bamCoverage` | Generate BigWig coverage tracks | `--binSize`, `--normalizeUsing`, `--ignoreDuplicates`, `-p` |
| `macs3` | Call peaks from aligned reads | `-t`, `-c`, `-f BAMPE`, `-g`, `-n`, `--outdir`, `-q` |
| `multiqc` | Aggregate QC reports | `-c` (config file) |
| `annotatePeaks.pl` | Annotate peaks to genomic features | `-gtf` (annotation file) |
| `computeMatrix` | Create matrix for heatmaps | `-S`, `-R`, `-a`, `-b`, `--referencePoint`, `--skipZeros`, `-p` | 
| `plotHeatmap` | Plot heatmaps from matrix | `-m`, `-o`, `--colorMap`, `--samplesLabel` |
---

### For Loops for Batch Processing

```bash
FILES=("file1" "file2" "file3")

for i in ${!FILES[@]}; do
    <commands using ${FILES[$i]}>
done
```

---

### Key Concepts Learned

**Paired-end sequencing:**
- Each DNA fragment sequenced from both ends (R1 and R2)
- Provides more accurate mapping and insert size information
- "Concordant pairs" = both mates align properly

**Quality control:**
- Performed automatically AFTER trimming via trim_galore
- FastQC reports show base quality, adapter content, duplication levels

**Normalization:**
- CPM (Counts Per Million) normalization allows fair comparison between samples
- Essential when samples have different sequencing depths

**Peak calling:**
- Treatment vs. control comparison identifies true enrichment
- Q-value cutoff controls false discovery rate
- H3K4me3 shows ~3000 peaks; DOT5-FLAG shows ~1 peak

**Biological interpretation:**
- H3K4me3 enriched at promoter-TSS (75% of peaks)
- Consistent with known role as active promoter mark
- DOT5 shows no specific genomic enrichment despite chromatin binding

---

### Cleanup and Shutdown

Check total disk usage:
```bash
du -sh ~/cutrun/
```

**Stop your EC2 instance** to avoid charges:
```bash
exit  # Exit SSH session
```

Then from AWS Console:
1. Go to EC2 Dashboard
2. Select your instance
3. Instance State → Stop instance

---

**Congratulations on completing the CUT&RUN analysis tutorial!**