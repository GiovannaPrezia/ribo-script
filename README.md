# RiboLongSmORF

Automated Ribo-seq processing and lncRNA-derived smORF discovery pipeline.

RiboLongSmORF provides an end-to-end workflow for processing public or in-house Ribo-seq datasets, from SRA download to alignment, quantification, quality control, and downstream analyses focused on small open reading frames (smORFs) encoded by long non-coding RNAs (lncRNAs).

---

## Features

* Automatic download of public datasets from NCBI SRA
* FASTQ generation using fasterq-dump
* FastQC and MultiQC reports
* Adapter trimming with Cutadapt
* Removal of contaminant RNAs (rRNA, tRNA, snRNA, snoRNA)
* STAR genome alignment
* featureCounts quantification
* Automated Ribo-seq quality control figures
* PCA analysis
* P-site regional distribution
* Metagene profiling
* Frame preference analysis
* 3-nt periodicity analysis
* Final QC reports and summary tables

---

## Workflow

```text
SRA Download
    ↓
FASTQ Generation
    ↓
FastQC (raw)
    ↓
Cutadapt Trimming
    ↓
FastQC (trimmed)
    ↓
noN / noPolyG Filtering
    ↓
Length Selection
    ├── all_lengths
    └── 28_36
            ↓
Contaminant Removal (Bowtie)
            ↓
FastQC (clean)
            ↓
STAR Alignment
            ↓
featureCounts
            ↓
MultiQC
            ↓
Ribo-seq QC Figures
            ↓
Final Reports
```

---

## Requirements

* Linux
* Conda (Miniconda or Anaconda)
* Internet connection for downloading reference files and SRA datasets

---

## Installation

Clone the repository:

```bash
git clone https://github.com/GiovannaPrezia/RiboLongSmorf.git

cd RiboLongSmorf
```

Run the setup and pipeline:

```bash
bash RiboLongSmORF.sh config.yaml
```

The script automatically:

* Creates the Conda environment
* Downloads reference files
* Builds the STAR genome index
* Starts the analysis pipeline

---

## Recommended: Use Screen

Large projects may run for several hours.

Create a screen session before starting:

```bash
screen -S ribolongsmorf
```

Run the pipeline:

```bash
bash RiboLongSmORF.sh config.yaml
```

Detach from the session:

```bash
Ctrl+A D
```

Reconnect later:

```bash
screen -r ribolongsmorf
```

List active sessions:

```bash
screen -ls
```

---

## Example Configuration

```yaml
project_name: California_PRJNA544411

project_description: Human cardiomyocyte differentiation Ribo-seq

project_root: /data/California_PRJNA544411

threads: 20

size_modes:
  - all_lengths
  - 28_36

samples:
  - run_id: SRR9113067
    sample_name: D15_rep1

  - run_id: SRR9113068
    sample_name: D15_rep2

  - run_id: SRR9113069
    sample_name: D15_rep3
```

---

## Output Structure

```text
01_SRAs/
02_fastq/
03_trimmed/
04_cleaned/
05_alignment/
06_star_qc/
07_counts/
08_annotation/
09_genome/
10_Ribotricer/
11_MultiQC/
12_QC_Figures/
13_Report/

logs/
QC_tables/
```

---

## Generated QC Figures

| Script | Description                  |
| ------ | ---------------------------- |
| 01     | Read length distribution     |
| 02     | P-site regional distribution |
| 03     | P-site metagene profile      |
| 04     | PCA                          |
| 05     | Frame preference             |
| 06     | 3-nt periodicity             |
| 07     | Alignment summary            |
| 08     | Contaminant summary          |

---

## References

Reference files are downloaded automatically:

* GRCh38 primary assembly genome
* GENCODE v45 annotation
* Human contaminant RNA database

---

Author: Giovanna N. B. Prezia
