#!/bin/bash

CONDA_ENV_PY3=encode-chip-seq-pipeline
CONDA_ENV_PY2=encode-chip-seq-pipeline-python2
CONDA_ENV_OLD_PY3=encode-chip-seq-pipeline-python3

conda env remove -n ${CONDA_ENV_PY3} -y
conda env remove -n ${CONDA_ENV_PY2} -y
conda env remove -n ${CONDA_ENV_OLD_PY3} -y

