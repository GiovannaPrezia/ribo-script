#!/bin/bash
set -e

CONFIG_FILE="${1:-config.yaml}"
ENV_NAME="ribolongsmorf_env"

echo "======================================"
echo "       RiboLongSmORF Pipeline"
echo "======================================"

# ======================================
# CONDA ENVIRONMENT
# ======================================

echo ""
echo "Checking Conda environment..."

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' found."
else
    echo "Environment '$ENV_NAME' not found."
    echo "Creating environment from environment.yml..."

    conda env create -f environment.yml

    echo ""
    echo "Environment successfully created."
fi

# ======================================
# REFERENCES
# ======================================

echo ""
echo "Downloading references..."

conda run -n "$ENV_NAME" \
    bash scripts/setup/01_download_references.sh "$CONFIG_FILE"

# ======================================
# STAR INDEX
# ======================================

echo ""
echo "Building STAR index..."

conda run -n "$ENV_NAME" \
    bash scripts/setup/02_build_star_index.sh "$CONFIG_FILE"

# ======================================
# PIPELINE
# ======================================

echo ""
echo "Starting pipeline..."

conda run -n "$ENV_NAME" \
    bash riboshortlong_pipe.sh "$CONFIG_FILE"
