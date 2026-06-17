# 🧬 RiboLongShort Pipeline: Detection of smORFs in lncRNAs from Ribo-seq Data

RiboLongShort is an automated pipeline designed for processing Ribo-seq datasets and identifying translated small open reading frames (smORFs), with a particular focus on lncRNA-derived smORFs.

The pipeline performs preprocessing, quality control, contaminant removal, genome alignment, quantification, and ORF detection, producing standardized outputs suitable for downstream translational analysis.

The final outputs include:

* Quality control reports
* Cleaned FASTQ files
* Genome-aligned BAM files
* Gene count matrices
* Ribotricer smORF predictions
* Integrated MultiQC reports

These outputs can be directly used for downstream analyses such as candidate prioritization, differential translation, functional annotation, visualization, and experimental validation.

---

## 🔁 Workflow Overview

The pipeline performs the following steps:

1. Download SRA data (optional)
2. Convert SRA files to FASTQ (`fasterq-dump`)
3. Raw quality control (`FastQC`)
4. Adapter trimming (`Cutadapt`)
5. Post-trimming quality assessment
6. Contaminant removal (`Bowtie1`)
7. Quality control of cleaned reads
8. Genome alignment (`STAR`)
9. Read quantification (`featureCounts`)
10. Alignment quality assessment
11. Integrated report generation (`MultiQC`)
12. smORF detection (`Ribotricer`)

---

## 📂 Project Structure

The repository provides an automatic project structure generator:

```bash
bash create_project_structure.sh
```

This creates the standardized directory layout required by the pipeline.

---

## ⚙️ Configuration

All project-specific settings are defined in:

```bash
config.yaml
```

The user only needs to edit:

* Project name
* Project location
* Resource directory
* Sample names
* SRA accessions
* Number of threads

No modifications to the pipeline code are required.

---

## 🧪 Software Requirements

The pipeline uses:

* SRA Toolkit
* FastQC
* Cutadapt
* Bowtie1
* STAR
* Samtools
* Subread (featureCounts)
* MultiQC
* Ribotricer

All dependencies are provided through the Docker image.

---

## 🐳 Docker

Build the container:

```bash
docker build -t ribolongshort .
```

Run the container:

```bash
docker run -it \
-v $(pwd):/workspace \
ribolongshort
```

---

## 🚀 Running the Pipeline

Edit:

```bash
config.yaml
```

Then run:

```bash
bash ribolongshort_pipe.sh config.yaml
```

The pipeline can be executed either:

* Step-by-step (interactive mode)
* End-to-end (full pipeline mode)

---

## 📊 Outputs

Main outputs include:

```text
02_fastq/
03_trimmed/
04_cleaned/
05_alignment/
06_star_qc/
07_counts/
10_Ribotricer/
11_MultiQC/
```

The final report contains:

* Sample information
* FastQC summaries
* Contaminant filtering statistics
* STAR alignment metrics
* featureCounts statistics
* Ribotricer summary statistics

---

## 🔬 Future Development

Planned additions include:

* iRibo integration
* Automated smORF prioritization
* ORF annotation modules
* R-based QC figure generation
* Snakemake implementation
* Comparative analysis modules
