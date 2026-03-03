#!/usr/bin/env bash
#
# Script: build_genome_index.sh
# Description:
#   Download GRCh38 genome (GENCODE), prepare STAR and Bowtie1 indices.
#
# Requirements:
#   - wget
#   - STAR
#   - bowtie (v1)
#
# Author: Giovanna N. B. Prezia
# License: MIT
#

set -euo pipefail

############################################
# USER CONFIGURATION
############################################

BASE_DIR="/path/to/project"
THREADS=50
READ_LENGTH=30              # Adjust to your read length
GENCODE_VERSION=45          # GENCODE release version
RAM_LIMIT=60000000000       # 60GB for STAR genomeGenerate

############################################
# DIRECTORIES
############################################

GENOME_DIR="${BASE_DIR}/09_genome"
ANNOTATION_DIR="${BASE_DIR}/08_annotation"
STAR_INDEX="${GENOME_DIR}/hg38_star_index"
BOWTIE1_INDEX="${GENOME_DIR}/hg38_bowtie1_index"

mkdir -p \
    "$GENOME_DIR" \
    "$ANNOTATION_DIR" \
    "$STAR_INDEX" \
    "$BOWTIE1_INDEX"

############################################
# CHECK DEPENDENCIES
############################################

command -v STAR >/dev/null 2>&1 || { echo "STAR not found."; exit 1; }
command -v bowtie-build >/dev/null 2>&1 || { echo "bowtie-build not found."; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "wget not found."; exit 1; }

############################################
# DOWNLOAD GENOME (GENCODE)
############################################

cd "$GENOME_DIR"

GENOME_FILE="GRCh38.primary_assembly.genome.fa"
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/${GENOME_FILE}.gz"

if [[ ! -f "$GENOME_FILE" ]]; then
    echo "Downloading GRCh38 genome (GENCODE v${GENCODE_VERSION})..."
    wget "$GENOME_URL"
    gunzip "${GENOME_FILE}.gz"
fi

GENOME_FA="${GENOME_DIR}/${GENOME_FILE}"
echo "Genome ready."

############################################
# DOWNLOAD ANNOTATION (GENCODE)
############################################

cd "$ANNOTATION_DIR"

GTF_FILE="gencode.v${GENCODE_VERSION}.annotation.gtf"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/${GTF_FILE}.gz"

if [[ ! -f "$GTF_FILE" ]]; then
    echo "Downloading GENCODE annotation v${GENCODE_VERSION}..."
    wget "$GTF_URL"
    gunzip "${GTF_FILE}.gz"
fi

GTF_PATH="${ANNOTATION_DIR}/${GTF_FILE}"
echo "Annotation ready."

############################################
# BUILD STAR INDEX
############################################

echo "Building STAR index..."

STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDEX} \
    --genomeFastaFiles ${GENOME_FA} \
    --sjdbGTFfile ${GTF_PATH} \
    --sjdbOverhang $((READ_LENGTH - 1)) \
    --limitGenomeGenerateRAM ${RAM_LIMIT}

echo "STAR index completed."

############################################
# BUILD BOWTIE1 INDEX
############################################

echo "Building Bowtie1 index..."

cd "${BOWTIE1_INDEX}"
bowtie-build ${GENOME_FA} hg38_index

echo "Bowtie1 index completed."

echo "All genome indices successfully generated."
