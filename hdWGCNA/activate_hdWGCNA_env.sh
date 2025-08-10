#!/bin/bash

set -e

# === CONFIG ===
MINICONDA_DIR="$HOME/miniconda3"
CONDA_ENV_NAME="hdwgcna_env"
PYTHON_VERSION="3.10"
R_VERSION="4.3"

echo "=== [1/7] Downloading Miniconda installer ==="
cd ~
wget -nc https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

echo "=== [2/7] Installing or Updating Miniconda at $MINICONDA_DIR ==="
if [ -d "$MINICONDA_DIR" ]; then
    echo "Miniconda already exists. Updating in place..."
    bash Miniconda3-latest-Linux-x86_64.sh -u -b -p $MINICONDA_DIR
else
    echo "Miniconda not found. Installing fresh..."
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINICONDA_DIR
fi

echo "=== [3/7] Initializing Conda ==="
eval "$($MINICONDA_DIR/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc

echo "=== [4/7] Accepting Conda Terms of Service ==="
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

echo "=== [5/7] Configuring Conda-forge channel ==="
conda config --add channels conda-forge
conda config --set channel_priority strict

echo "=== [6/7] Creating Conda environment: $CONDA_ENV_NAME ==="
conda create -y -n $CONDA_ENV_NAME python=$PYTHON_VERSION r-base=$R_VERSION r-essentials

echo "=== [7/7] Installing Python and R dependencies ==="
conda activate $CONDA_ENV_NAME

# Python packages
conda install -y numpy pandas scipy scikit-learn matplotlib seaborn umap-learn anndata scanpy

# Core R packages from conda-forge
conda install -y r-reticulate r-devtools r-tidyverse r-data.table r-remotes

echo "=== [8/8] Installing R packages from GitHub (specific versions) ==="
# UCell (required for hdWGCNA)
Rscript -e "remotes::install_github('carmonalab/UCell')"

# Seurat v4.4.0 (!!!) specifically
Rscript -e "remotes::install_github('satijalab/seurat', ref = 'v4.4.0')"

# hdWGCNA
Rscript -e "devtools::install_github('smorabit/hdWGCNA', ref = 'dev')"

# RedRibbon
Rscript -e "devtools::install_github('antpiron/RedRibbon')"


echo "All done!"
echo "To activate your environment in future sessions, run:"
echo "    conda activate $CONDA_ENV_NAME"
