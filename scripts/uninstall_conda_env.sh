#!/bin/bash

PIPELINE_CONDA_ENVS=(
  encd-chip
  encd-chip-macs2
  encd-chip-spp
)
for PIPELINE_CONDA_ENV in "${PIPELINE_CONDA_ENVS[@]}"
do
  conda env remove -n ${PIPELINE_CONDA_ENV} -y
done
