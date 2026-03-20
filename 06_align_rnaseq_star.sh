#!/bin/bash
# STAR alignment — RNA-seq (paired-end)

set -euo pipefail

THREADS=20

ulimit -n 10000

############################################
# DIRECTORIES
############################################

BASE_DIR="/home/giovanna.prezia/Diretório/data/bioprojects/Columbia_Data"

RNA_DIR="${BASE_DIR}/02_fastq/rna_seq"
ALIGN_DIR="${BASE_DIR}/05_alignment/rna_seq"
STAR_QC="${BASE_DIR}/06_star_qc"

STAR_INDEX="${BASE_DIR}/09_genome/hg38_star_index"

############################################
# RNASEQ ALIGNMENT
############################################

echo "=============================="
echo "RNA-seq STAR alignment"
echo "=============================="

shopt -s nullglob

for R1 in "${RNA_DIR}"/*_1.fastq.gz
do

    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    SAMPLE=$(basename "${R1}" _1.fastq.gz)

    if [ ! -f "${R2}" ]; then
        echo "⚠️ Missing pair for ${SAMPLE}"
        continue
    fi

    BAM="${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

    if [ -f "${BAM}" ]; then
        echo "⏩ ${SAMPLE} already aligned"
        continue
    fi

    echo "🚀 Aligning ${SAMPLE}"

    STAR \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMstrandField intronMotif \
        --outFilterMultimapNmax 10


    samtools flagstat -@ ${THREADS} "${BAM}" \
        > "${STAR_QC}/${SAMPLE}_rnaseq_flagstat.txt"

    echo "✅ Done: ${SAMPLE}"
    echo

done

echo "🎉 RNA-seq alignment finished"
