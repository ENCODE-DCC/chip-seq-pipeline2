#!/bin/bash

PIPELINE_CONDA_ENVS=(
  encode-chip-seq-pipeline
  encode-chip-seq-pipeline-macs2
  encode-chip-seq-pipeline-spp
)
for PIPELINE_CONDA_ENV in "${PIPELINE_CONDA_ENVS[@]}"
do
  conda env remove -n ${PIPELINE_CONDA_ENV} -y
done
