#!/bin/bash

# Paths and filenames
LOCAL_DIR="$HOME/km_MousePOASnseq/Scripts/km_cleanScripts/data"
CONTAINER_DIR="/data"
SCRIPT_DIR="/km_rscvi_docker"
SCRIPT_NAME="run_mapscvi.R"

# Run the container and execute the R script
docker run --rm \
  -v "${LOCAL_DIR}:${CONTAINER_DIR}" \
  km_rmapscvi \
  Rscript "${CONTAINER_DIR}${SCRIPT_DIR}/${SCRIPT_NAME}"
