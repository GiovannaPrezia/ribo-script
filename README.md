# 🧬 RiboLongShort Pipeline: Detection of smORFs in lncRNAs from riboseq data 

This repository contains an automated **Bash-based pipeline** for processing **Ribo-seq** data, with a focus on detecting **small open reading frames (smORFs)** encoded by **long non-coding RNAs (lncRNAs)**.  
The pipeline includes full preprocessing of raw sequencing data and identification of translated ORFs using **Ribotricer**.

The final output includes aligned BAM files, count matrices, QC reports, and Ribotricer-generated smORF predictions ready for downstream analysis and interpretation.

---

## 🔁 Workflow Overview

The pipeline runs the following steps:

1. **Download SRA files** using `prefetch`  
2. **Convert SRA to FASTQ** using `fasterq-dump`  
3. **Quality control** using **FastQC**  
4. **Adapter trimming** with **Cutadapt**  
5. **Contaminant removal** (rRNA, tRNA, etc.) using **Bowtie1**  
6. **Genome alignment** using **STAR**  
7. **Read quantification** using **featureCounts**  
8. **Alignment QC** using STAR log files  
9. **Integrated quality report** using **MultiQC**  
10. **smORF detection** using **Ribotricer** (frame periodicity + ribosome occupancy)

All steps are orchestrated by the master script `pipeline_ribo-seq.sh`, with modular Bash scripts stored in the `scripts/` directory.

---

## 🧪 Requirements

The following tools must be installed (preferably in a dedicated conda environment):

- `sra-tools`
- `fastqc`
- `cutadapt`
- `bowtie` (Bowtie1)
- `STAR`
- `samtools`
- `subread` (featureCounts)
- `multiqc`
- `ribotricer` (Python package)

