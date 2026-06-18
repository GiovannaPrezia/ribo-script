# 🧬 RiboLongSmORF Pipeline: Detection of smORFs in lncRNAs from Ribo-seq Data

RiboLongSmORF is an automated workflow for processing Ribo-seq datasets, performing comprehensive quality control, and identifying translated small open reading frames (smORFs), with a particular focus on lncRNA-derived smORFs and putative microproteins.

The pipeline integrates preprocessing, quality assessment, contaminant removal, genome alignment, read quantification, and ORF detection into a reproducible and standardized workflow suitable for translational profiling studies.

The final outputs include:

* Quality control reports
* Cleaned FASTQ files
* Genome-aligned BAM files
* Gene count matrices
* Ribotricer smORF predictions
* Integrated MultiQC reports
* Ribo-seq-specific QC figures

These outputs can be directly used for downstream analyses such as candidate prioritization, differential translation studies, functional annotation, visualization, and experimental validation.

---

## 🔁 Workflow Overview

RiboLongSmORF performs the following steps:

1. Download SRA data (optional)
2. Convert SRA files to FASTQ (`fasterq-dump`)
3. Raw read quality control (`FastQC`)
4. Adapter and sequence trimming (`Cutadapt`)
5. Post-trimming quality assessment
6. Contaminant removal (`Bowtie1`)
7. Quality control of cleaned reads
8. Genome alignment (`STAR`)
9. Read quantification (`featureCounts`)
10. Alignment quality assessment
11. Integrated report generation (`MultiQC`)
12. Ribo-seq QC figure generation
13. smORF detection (`Ribotricer`)

---

## ✨ Key Features

* Automated end-to-end Ribo-seq processing
* Support for public SRA datasets
* Contaminant filtering using a curated human RNA contaminant database
* Standardized STAR alignment and featureCounts quantification
* Integrated MultiQC reports
* Ribo-seq-specific quality control metrics
* Automated smORF identification with Ribotricer
* Optimized for lncRNA-derived smORF and microprotein discovery
* Reproducible and configurable workflow through a single `config.yaml` file
