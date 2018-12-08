# ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline

[![CircleCI](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master)

## Introduction 
This ChIP-Seq pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications (by Anshul Kundaje) in [this google doc](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#).

### Features

* **Flexibility**: Support for `docker`, `singularity` and `Conda`.
* **Portability**: Support for many cloud platforms (Google/DNANexus) and cluster engines (SLURM/SGE/PBS).
* **Beatiful HTML report**: tabulated quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures).
  - Examples: [HTML](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/example_output/qc.html), [JSON](docs/example_output/qc.json)
* **Genomes**: Pre-built database for GRCh38, hg19, mm10, mm9 and additional support for custom genomes.

## Installation and tutorial

This pipeline supports many cloud platforms and cluster engines. It also supports `docker`, `singularity` and `Conda` to resolve complicated software dependencies for the pipeline. A tutorial-based instruction for each platform will be helpful to understand how to run pipelines. There are special instructions for two major Stanford HPC servers (SCG4 and Sherlock).

* Cloud platforms
  * Web interface
    * [DNANexus Platform](docs/tutorial_dx_web.md)
  * CLI (command line interface)
    * [Google Cloud Platform](docs/tutorial_google.md)
    * [DNANexus Platform](docs/tutorial_dx_cli.md)
* Stanford HPC servers (CLI)
  * [Stanford SCG4](docs/tutorial_scg.md)
  * [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)
* Cluster engines (CLI)
  * [SLURM](docs/tutorial_slurm.md)
  * [Sun GridEngine (SGE/PBS)](docs/tutorial_sge.md)
* Local computers (CLI)
  * [Local system with `singularity`](docs/tutorial_local_singularity.md)
  * [Local system with `docker`](docs/tutorial_local_docker.md)
  * [Local system with `Conda`](docs/tutorial_local_conda.md)

## Input JSON file

[Input JSON file specification](docs/input.md)

## Output directories

[Output directory specification](docs/output.md)