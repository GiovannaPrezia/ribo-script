#!/usr/bin/env bash
#
# Filter reads mapping to abundant structural / non-target RNAs (e.g., rRNA, tRNA, snoRNA)
# from trimmed Ribo-seq reads using Bowtie1.
#
# The script keeps reads that do NOT align to this reference set (via --un).
#
# Author: Giovanna N. B. Prezia
#

set -euo pipefail

# -----------------------------
# USER PARAMETERS (EDIT HERE)
# -----------------------------
BASE_DIR="/home/giovanna.prezia/Diretório/data/bioprojects/Columbia_Data"
THREADS=20
MISMATCHES=0

TRIMMED_DIR="${BASE_DIR}/03_trimmed/ribo_seq"
CLEANED_DIR="${BASE_DIR}/04_cleaned/ribo_seq"

# Bowtie1 index prefix for your structural/non-target RNAs set
# Must point to files like: <prefix>.1.ebwt, <prefix>.2.ebwt, etc.
INDEX_PREFIX="${BASE_DIR}/09_genome/rnas_dictionary/bowtie1_index_clean/rna_dictionary_bt1"

# -----------------------------
# CHECKS
# -----------------------------
command -v bowtie >/dev/null 2>&1 || { echo "ERROR: bowtie not found in PATH"; exit 1; }
command -v pigz  >/dev/null 2>&1 || { echo "ERROR: pigz not found in PATH"; exit 1; }

if [[ ! -f "${INDEX_PREFIX}.1.ebwt" ]]; then
  echo "ERROR: Bowtie1 index not found: ${INDEX_PREFIX}.1.ebwt"
  echo "Check INDEX_PREFIX."
  exit 1
fi

mkdir -p "${CLEANED_DIR}"

shopt -s nullglob
INPUTS=( "${TRIMMED_DIR}"/*.trim.fastq.gz "${TRIMMED_DIR}"/*.trim.fastq )

if [[ ${#INPUTS[@]} -eq 0 ]]; then
  echo "ERROR: No input files found in ${TRIMMED_DIR}"
  echo "Expected: *.trim.fastq.gz or *.trim.fastq"
  exit 1
fi

echo "Filtering reads mapping to structural / non-target RNAs with Bowtie1"
echo "Index:  ${INDEX_PREFIX}"
echo "Input:  ${TRIMMED_DIR}"
echo "Output: ${CLEANED_DIR}"
echo "Threads: ${THREADS}"
echo

# -----------------------------
# FUNCTIONS
# -----------------------------
run_bowtie_gzip() {
  local infile="$1"
  local out_un="$2"
  local log="$3"

  bowtie \
    --gzip \
    -q \
    -p "${THREADS}" \
    -v "${MISMATCHES}" \
    --un "${out_un}" \
    "${INDEX_PREFIX}" \
    "${infile}" \
    > /dev/null 2> "${log}"
}

run_bowtie_stdin() {
  local infile="$1"
  local out_un="$2"
  local log="$3"

  # Use pigz to stream-decompress; more robust than zcat on some systems
  pigz -dc "${infile}" | bowtie \
    -q \
    -p "${THREADS}" \
    -v "${MISMATCHES}" \
    --un "${out_un}" \
    "${INDEX_PREFIX}" \
    - \
    > /dev/null 2> "${log}"
}

run_bowtie_tempfastq() {
  local infile="$1"
  local out_un="$2"
  local log="$3"

  local tmpdir
  tmpdir="$(mktemp -d)"
  local tmpfastq="${tmpdir}/input.fastq"

  pigz -dc "${infile}" > "${tmpfastq}"

  bowtie \
    -q \
    -p "${THREADS}" \
    -v "${MISMATCHES}" \
    --un "${out_un}" \
    "${INDEX_PREFIX}" \
    "${tmpfastq}" \
    > /dev/null 2> "${log}"

  rm -rf "${tmpdir}"
}

# -----------------------------
# MAIN LOOP
# -----------------------------
for FILE in "${INPUTS[@]}"; do
  base="$(basename "${FILE}")"
  sample="${base%.trim.fastq.gz}"
  sample="${sample%.trim.fastq}"

  OUT_UN="${CLEANED_DIR}/${sample}.clean.fastq"
  OUT_UN_GZ="${OUT_UN}.gz"
  LOG="${CLEANED_DIR}/${sample}.bowtie1.log"

  if [[ -f "${OUT_UN_GZ}" ]]; then
    echo "✔ ${sample} already processed. Skipping."
    continue
  fi

  echo "Processing ${sample}..."
  echo "  Input: ${FILE}"

  # Remove any partial outputs from previous failed runs
  rm -f "${OUT_UN}" "${OUT_UN_GZ}" "${LOG}"

  status=0

  if [[ "${FILE}" == *.gz ]]; then
    # Attempt 1: bowtie --gzip
    set +e
    run_bowtie_gzip "${FILE}" "${OUT_UN}" "${LOG}"
    status=$?
    set -e

    if [[ $status -ne 0 ]]; then
      echo "  ⚠ Attempt 1 failed (bowtie --gzip). Trying stdin streaming..."

      # Attempt 2: streaming stdin
      rm -f "${OUT_UN}"
      set +e
      run_bowtie_stdin "${FILE}" "${OUT_UN}" "${LOG}"
      status=$?
      set -e
    fi

    if [[ $status -ne 0 ]]; then
      echo "  ⚠ Attempt 2 failed (stdin). Trying temporary FASTQ (guaranteed)..."

      # Attempt 3: temp fastq
      rm -f "${OUT_UN}"
      set +e
      run_bowtie_tempfastq "${FILE}" "${OUT_UN}" "${LOG}"
      status=$?
      set -e
    fi

  else
    # Non-gz input: run directly
    set +e
    bowtie \
      -q \
      -p "${THREADS}" \
      -v "${MISMATCHES}" \
      --un "${OUT_UN}" \
      "${INDEX_PREFIX}" \
      "${FILE}" \
      > /dev/null 2> "${LOG}"
    status=$?
    set -e
  fi

  if [[ $status -ne 0 ]]; then
    echo "❌ Bowtie failed for ${sample} (exit code: ${status})"
    echo "---- log (last 60 lines): ${LOG} ----"
    tail -n 60 "${LOG}" || true
    echo "-------------------------------------"
    exit 1
  fi

  if [[ ! -s "${OUT_UN}" ]]; then
    echo "❌ Output file was not created or is empty: ${OUT_UN}"
    echo "---- log (last 60 lines): ${LOG} ----"
    tail -n 60 "${LOG}" || true
    echo "-------------------------------------"
    exit 1
  fi

  pigz -f -p "${THREADS}" "${OUT_UN}"
  echo "✅ Clean reads: ${OUT_UN_GZ}"
  echo "📄 Log: ${LOG}"
  echo
done

echo "🎉 Done."
