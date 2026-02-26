#!/usr/bin/env bash
#
# Script: download_srr_fastq_fastqc.sh
# Description:
#   Download SRA accessions, convert to FASTQ, auto-fix unexpected .fast extensions, compress with pigz and run FastQC.
#
# Requirements:
#   - SRA Toolkit (prefetch, fasterq-dump)
#   - FastQC
#   - pigz
#
# Author: Your Name
# License: MIT
#

set -euo pipefail

############################################
# USAGE
############################################

usage() {
    echo "Usage:"
    echo "  $0 -b BASE_DIR -r RIBO_LIST -n RNA_LIST -t THREADS"
    echo ""
    echo "Arguments:"
    echo "  -b   Base project directory"
    echo "  -r   File with Ribo-seq SRR accessions (single-end)"
    echo "  -n   File with RNA-seq SRR accessions (paired-end)"
    echo "  -t   Number of threads (default: 8)"
    echo ""
    exit 1
}

############################################
# DEFAULTS
############################################

THREADS=8

############################################
# ARGUMENT PARSING
############################################

while getopts "b:r:n:t:h" opt; do
    case ${opt} in
        b) BASE_DIR=$OPTARG ;;
        r) RIBO_SRR_LIST=$OPTARG ;;
        n) RNA_SRR_LIST=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${BASE_DIR:-}" || -z "${RIBO_SRR_LIST:-}" || -z "${RNA_SRR_LIST:-}" ]]; then
    usage
fi

############################################
# DIRECTORIES
############################################

SRA_DIR="${BASE_DIR}/01_SRAs"
FASTQ_DIR="${BASE_DIR}/02_fastq"

mkdir -p \
    "${SRA_DIR}/ribo_seq" \
    "${SRA_DIR}/rna_seq" \
    "${FASTQ_DIR}/ribo_seq/fastqc_raw" \
    "${FASTQ_DIR}/rna_seq/fastqc_raw"

############################################
# FUNCTION: Auto-fix .fast extension
############################################

fix_fast_extension() {
    local FILE_BASE=$1

    if [[ -f "${FILE_BASE}.fast" && ! -f "${FILE_BASE}.fastq" ]]; then
        echo "⚠ Fixing extension for ${FILE_BASE}"
        mv "${FILE_BASE}.fast" "${FILE_BASE}.fastq"
    fi
}

############################################
# RIBO-SEQ (Single-End)
############################################

echo "=== Starting Ribo-seq processing ==="

while read -r SRR || [[ -n "$SRR" ]]; do
    [[ -z "$SRR" ]] && continue

    FINAL_FASTQ="${FASTQ_DIR}/ribo_seq/${SRR}.fastq.gz"

    if [[ -f "${FINAL_FASTQ}" ]]; then
        echo "✔ ${SRR} already processed"
        continue
    fi

    if [[ ! -f "${SRA_DIR}/ribo_seq/${SRR}/${SRR}.sra" ]]; then
        prefetch "${SRR}" --output-directory "${SRA_DIR}/ribo_seq"
    fi

    fasterq-dump \
        "${SRA_DIR}/ribo_seq/${SRR}/${SRR}.sra" \
        -O "${FASTQ_DIR}/ribo_seq" \
        --threads "${THREADS}"

    fix_fast_extension "${FASTQ_DIR}/ribo_seq/${SRR}"

    pigz -f -p "${THREADS}" "${FASTQ_DIR}/ribo_seq/${SRR}.fastq"

    fastqc "${FINAL_FASTQ}" \
        -o "${FASTQ_DIR}/ribo_seq/fastqc_raw" \
        --threads "${THREADS}"

done < "${RIBO_SRR_LIST}"

############################################
# RNA-SEQ (Paired-End)
############################################

echo "=== Starting RNA-seq processing ==="

while read -r SRR || [[ -n "$SRR" ]]; do
    [[ -z "$SRR" ]] && continue

    R1="${FASTQ_DIR}/rna_seq/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/rna_seq/${SRR}_2.fastq.gz"

    if [[ -f "${R1}" && -f "${R2}" ]]; then
        echo "✔ ${SRR} already processed"
        continue
    fi

    if [[ ! -f "${SRA_DIR}/rna_seq/${SRR}/${SRR}.sra" ]]; then
        prefetch "${SRR}" --output-directory "${SRA_DIR}/rna_seq"
    fi

    fasterq-dump \
        "${SRA_DIR}/rna_seq/${SRR}/${SRR}.sra" \
        -O "${FASTQ_DIR}/rna_seq" \
        --split-files \
        --threads "${THREADS}"

    for read in 1 2; do
        fix_fast_extension "${FASTQ_DIR}/rna_seq/${SRR}_${read}"
    done

    pigz -f -p "${THREADS}" "${FASTQ_DIR}/rna_seq/${SRR}_1.fastq"
    pigz -f -p "${THREADS}" "${FASTQ_DIR}/rna_seq/${SRR}_2.fastq"

    fastqc "${R1}" "${R2}" \
        -o "${FASTQ_DIR}/rna_seq/fastqc_raw" \
        --threads "${THREADS}"

done < "${RNA_SRR_LIST}"

echo "=== Pipeline completed successfully ==="
