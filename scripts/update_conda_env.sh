#!/bin/bash
set -e  # Stop on error

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
SRC_DIR=${SH_SCRIPT_DIR}/../src

PIPELINE_CONDA_ENVS=(
  encode-chip-seq-pipeline
  encode-chip-seq-pipeline-macs2
  encode-chip-seq-pipeline-spp
)
chmod u+rx ${SRC_DIR}/*.py

echo "** Updating WDL task wrappers on each Conda environment..."
for PIPELINE_CONDA_ENV in "${PIPELINE_CONDA_ENVS[@]}"
do	
  CONDA_BIN=$(dirname $(conda run -n ${PIPELINE_CONDA_ENV} which python))
  echo -e "** Transferring WDL task wrappers to ${CONDA_BIN}..."
  cp -f ${SRC_DIR}/*.py ${CONDA_BIN}/
done
