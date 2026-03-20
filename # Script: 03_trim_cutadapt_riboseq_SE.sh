#!/usr/bin/env bash
#
# Script: 03_trim_cutadapt_riboseq_SE.sh
#
# Purpose:
#   Pre-process Ribo-seq (single-end) reads using cutadapt
#   and generate FastQC reports after trimming.
#
# Input:
#   *.fastq.gz or *.fastq in IN_RIBO
#
# Output:
#   trimmed FASTQ
#   cutadapt log
#   FastQC report
#
# Author: Giovanna N. B. Prezia
#

set -euo pipefail

############################################
# USER SETTINGS
############################################

THREADS=20

BASE_DIR="/home/giovanna.prezia/Diretório/data/bioprojects/Columbia_Data"

IN_RIBO="${BASE_DIR}/02_fastq/ribo_seq"
OUT_RIBO="${BASE_DIR}/03_trimmed/ribo_seq"

FASTQC_DIR="${OUT_RIBO}/fastqc_trimmed"

mkdir -p "${OUT_RIBO}"
mkdir -p "${FASTQC_DIR}"

############################################
# DEPENDENCY CHECK
############################################

command -v cutadapt >/dev/null || { echo "ERROR: cutadapt not found."; exit 1; }
command -v fastqc >/dev/null || { echo "ERROR: fastqc not found."; exit 1; }

############################################
# INPUT FILES
############################################

shopt -s nullglob

FILES=( "${IN_RIBO}"/*.fastq.gz "${IN_RIBO}"/*.fastq )

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "ERROR: No FASTQ files found in ${IN_RIBO}"
  exit 1
fi

echo "======================================"
echo "Ribo-seq trimming with Cutadapt"
echo "======================================"

echo "Input:  ${IN_RIBO}"
echo "Output: ${OUT_RIBO}"
echo "FastQC: ${FASTQC_DIR}"
echo "Threads: ${THREADS}"
echo

############################################
# MAIN LOOP
############################################

for IN_FQ in "${FILES[@]}"; do

  base="$(basename "$IN_FQ")"
  sample="${base%.fastq.gz}"
  sample="${sample%.fastq}"

  OUT_FQ="${OUT_RIBO}/${sample}.trim.fastq.gz"
  LOG="${OUT_RIBO}/${sample}.cutadapt.log"

  FASTQC_HTML="${FASTQC_DIR}/${sample}.trim_fastqc.html"

  if [[ -f "${OUT_FQ}" ]]; then
    echo "✔ ${sample} already trimmed."
  else

    echo "Trimming ${sample}..."

    cutadapt \
      --cut 3 \
      -a "A{15}" \
      --nextseq-trim=20 \
      --minimum-length 17 \
      -j "${THREADS}" \
      -o "${OUT_FQ}" \
      "${IN_FQ}" \
      > "${LOG}" 2>&1

    echo "✅ Trimmed: ${OUT_FQ}"

  fi

  ########################################
  # FASTQC AFTER TRIMMING
  ########################################

  if [[ -f "${FASTQC_HTML}" ]]; then
    echo "✔ FastQC already exists for ${sample}"
  else

    echo "Running FastQC for ${sample}..."

    fastqc \
      -t "${THREADS}" \
      -o "${FASTQC_DIR}" \
      "${OUT_FQ}"

  fi

  echo

done

echo "======================================"
echo "Ribo-seq trimming completed"
echo "======================================"
