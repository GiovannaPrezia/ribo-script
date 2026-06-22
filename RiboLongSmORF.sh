#!/bin/bash
set -e

# ==========================================================
# CHECK CONDA
# ==========================================================

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found in PATH."
    echo "Please install Miniconda or Anaconda first."
    exit 1
fi

CONFIG_FILE="${1:-config.yaml}"
ENV_NAME="ribolongsmorf_env"

# ==========================================================
# CHECK FILES
# ==========================================================

[[ -f "$CONFIG_FILE" ]] || {
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
}

[[ -f "environment.yml" ]] || {
    echo "ERROR: environment.yml not found."
    exit 1
}

[[ -f "ribolongsmorf_pipe.sh" ]] || {
    echo "ERROR: ribolongsmorf_pipe.sh not found."
    exit 1
}

[[ -f "scripts/setup/01_download_references.sh" ]] || {
    echo "ERROR: 01_download_references.sh not found."
    exit 1
}

[[ -f "scripts/setup/02_build_star_index.sh" ]] || {
    echo "ERROR: 02_build_star_index.sh not found."
    exit 1
}

# ==========================================================
# HEADER
# ==========================================================

cat << EOF

╔════════════════════════════════════════════════════════════╗
║                      RiboLongSmORF                        ║
║                                                          ║
║     Ribo-seq Processing & lncRNA-smORF Discovery         ║
║                                                          ║
║     Automated pipeline for translational profiling       ║
╚════════════════════════════════════════════════════════════╝

EOF

# ==========================================================
# CONDA ENVIRONMENT
# ==========================================================

echo ""
echo "[1/4] Checking Conda environment..."
echo ""

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then

    echo "✓ Environment '$ENV_NAME' found."

else

    echo "Environment '$ENV_NAME' not found."
    echo "Creating environment from environment.yml..."
    echo ""

    conda env create -f environment.yml

    echo ""
    echo "✓ Environment successfully created."

fi

# ==========================================================
# REFERENCES
# ==========================================================

echo ""
echo "[2/4] Downloading reference files..."
echo ""

conda run -n "$ENV_NAME" \
    bash scripts/setup/01_download_references.sh "$CONFIG_FILE"

echo ""
echo "✓ Reference setup completed."

# ==========================================================
# STAR INDEX
# ==========================================================

echo ""
echo "[3/4] Building STAR genome index..."
echo ""

conda run -n "$ENV_NAME" \
    bash scripts/setup/02_build_star_index.sh "$CONFIG_FILE"

echo ""
echo "✓ STAR index ready."

# ==========================================================
# PIPELINE
# ==========================================================

echo ""
echo "[4/4] Starting RiboLongSmORF pipeline..."
echo ""

conda run -n "$ENV_NAME" \
    bash ribolongsmorf_pipe.sh "$CONFIG_FILE"

echo ""
echo "✓ RiboLongSmORF execution finished."
echo ""
