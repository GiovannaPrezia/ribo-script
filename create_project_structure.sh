#!/bin/bash
# create_project_rb_structure.sh
# Autor: Giovanna N. B. Prezia

set -e

echo "Criando estrutura organizada Ribo-seq only"

mkdir -p 01_SRAs/ribo_seq
mkdir -p 02_fastq/ribo_seq/fastqc_raw
mkdir -p 03_trimmed/ribo_seq/fastqc_trimmed
mkdir -p 04_cleaned/ribo_seq/fastqc_cleaned
mkdir -p 05_alignment/ribo_seq
mkdir -p 06_star_qc/ribo_seq
mkdir -p 07_counts/ribo_seq
mkdir -p 10_Ribotricer
mkdir -p 11_MultiQC
mkdir -p 12_scripts
mkdir -p logs
mkdir -p QC_tables

echo "Estrutura criada com sucesso 🧬"
