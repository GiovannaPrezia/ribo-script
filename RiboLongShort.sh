#!/bin/bash
set -e

CONFIG_FILE="${1:-config.yaml}"

echo "======================================"
echo "        RiboLongShort Pipeline"
echo "======================================"

echo "Checking conda environment..."


echo "Downloading references..."
bash scripts/setup/01_download_reference.sh "$CONFIG_FILE"

echo "Building STAR index..."
bash scripts/setup/02_build_star_index.sh "$CONFIG_FILE"

echo "Starting pipeline..."
bash ribolongshort_pipe.sh "$CONFIG_FILE"
