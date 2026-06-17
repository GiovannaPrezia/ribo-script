# 🧬 RiboLongShort Pipeline: Detection of smORFs in lncRNAs from Ribo-seq data

RiboLongShort is a modular pipeline for processing Ribo-seq datasets and identifying translated small open reading frames (smORFs), including candidate smORFs encoded by long non-coding RNAs (lncRNAs).

The workflow performs quality control, contaminant removal, genome alignment, quantification, and ORF detection using Ribotricer.

Features
SRA → FASTQ conversion
FastQC quality control
Adapter trimming with Cutadapt
rRNA/tRNA contaminant removal using Bowtie1
Genome alignment using STAR
Quantification using featureCounts
Integrated MultiQC reports
ORF detection using Ribotricer
Interactive or full-workflow execution
Workflow
SRA
 └── FASTQ
      └── FastQC
           └── Cutadapt
                └── Contaminant Removal (Bowtie1)
                     └── STAR Alignment
                          └── featureCounts
                               └── MultiQC
                                    └── Ribotricer
Requirements
SRA Toolkit
FastQC
MultiQC
Cutadapt
Bowtie1
STAR
samtools
Subread (featureCounts)
Ribotricer
Installation
git clone https://github.com/USERNAME/RiboLongShort.git
cd RiboLongShort
Usage

Run the interactive workflow:

bash scripts/run_riboseq_interactive.sh

The pipeline can be executed:

step-by-step (interactive mode)
fully automated (continuous mode)
Outputs

The pipeline generates:

Trimmed FASTQ files
Contaminant-filtered FASTQ files
Aligned BAM files
Alignment statistics
featureCounts tables
FastQC reports
MultiQC reports
Ribotricer ORF predictions
Roadmap
Configuration files (YAML)
Automated Ribotricer integration
QC figure generation in R
Snakemake workflow
Docker/Singularity support
iRibo integration


