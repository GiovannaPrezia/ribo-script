#!/bin/bash
# STAR alignment — Ribo-seq (single-end)

set -euo pipefail

THREADS=20

ulimit -n 10000

############################################
# DIRECTORIES
############################################

BASE_DIR="/home/giovanna.prezia/Diretório/data/bioprojects/Columbia_Data"

RIBO_DIR="${BASE_DIR}/04_cleaned/ribo_seq"
ALIGN_DIR="${BASE_DIR}/05_alignment/ribo_seq"
STAR_QC="${BASE_DIR}/06_star_qc"

STAR_INDEX="${BASE_DIR}/09_genome/hg38_star_index"

mkdir -p "${ALIGN_DIR}"
mkdir -p "${STAR_QC}"

############################################
# RIBOSEQ ALIGNMENT
############################################

echo "=============================="
echo "Ribo-seq STAR alignment"
echo "=============================="

shopt -s nullglob

for FILE in "${RIBO_DIR}"/*.clean.fastq.gz
do

    SAMPLE=$(basename "${FILE}" .clean.fastq.gz)

    BAM="${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

    if [ -f "${BAM}" ]; then
        echo "⏩ ${SAMPLE} already aligned"
        continue
    fi

    echo "🚀 Aligning ${SAMPLE}"

    STAR \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${FILE}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --alignSJoverhangMin 400 \
        --outFilterMismatchNmax 0 \
        --outFilterMatchNmin 15 \
        --outFilterMultimapNmax 1 \
        --quantMode TranscriptomeSAM

    samtools index -@ ${THREADS} "${BAM}"

    samtools flagstat -@ ${THREADS} "${BAM}" \
        > "${STAR_QC}/${SAMPLE}_riboseq_flagstat.txt"

    echo "✅ Done: ${SAMPLE}"
    echo

done

echo "🎉 Ribo-seq alignment finished"
