#!/bin/bash
set -e

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found in PATH."
    echo "Please install Miniconda or Anaconda first."
    exit 1
fi

CONFIG_FILE="${1:-config.yaml}"
ENV_NAME="ribolongsmorf_env"

[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config file not found: $CONFIG_FILE"; exit 1; }
[[ -f "environment.yml" ]] || { echo "ERROR: environment.yml not found."; exit 1; }
[[ -f "ribolongsmorf_pipe.sh" ]] || { echo "ERROR: ribolongsmorf_pipe.sh not found."; exit 1; }

[[ -f "scripts/setup/01_download_references.sh" ]] || { echo "ERROR: 01_download_references.sh not found."; exit 1; }
[[ -f "scripts/setup/02_build_star_index.sh" ]] || { echo "ERROR: 02_build_star_index.sh not found."; exit 1; }


cat << "EOF"

============================================================
                    RiboLongSmORF
------------------------------------------------------------
     Ribo-seq Processing and lncRNA-smORF Discovery
============================================================

EOF

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
    bash ribolongsmorf_pipe.sh "$CONFIG_FILE"
