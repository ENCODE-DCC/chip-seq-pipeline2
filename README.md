# ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline

[![CircleCI](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master)

## Introduction 
This ChIP-Seq pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications (by Anshul Kundaje) in [this google doc](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#).

### Features

* **Portability**: The pipeline run can be performed across different cloud platforms such as Google, AWS and DNAnexus, as well as on cluster engines such as SLURM, SGE and PBS.
* **User-friendly HTML report**: tabulated quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures).
  - Examples: [HTML](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/example_output/qc.html), [JSON](docs/example_output/v1.1.5/qc.json)
* **Supported genomes**: Pipeline needs genome specific data such as aligner indices, chromosome sizes file and blacklist. We provide a genome database downloader/builder for hg38, hg19, mm10, mm9. You can also use this [builder](docs/build_genome_database.md) to build genome database from FASTA for your custom genome.

## Installation
1) Git clone this repo.

	```bash
	$ cd
	$ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
	```

2) Install [Caper](https://github.com/ENCODE-DCC/caper#installation). Caper is a python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell).

	> **IMPORTANT**: Make sure that you have python3(> 3.4.1) installed on your system.

	```bash
	$ pip install caper  # use pip3 if it doesn't work
	```

3) Read through [Caper's README](https://github.com/ENCODE-DCC/caper) carefully. Find an instruction for your platform. 
	> **IMPORTANT**: Configure your Caper configuration file `~/.caper/default.conf` correctly for your platform.

## Running a pipeline locally with Caper

1) Prepare an input JSON file. We will use a subsampled example input JSON based on URLs. Caper will automatically download all fastqs and reference human genome data recursively.
	```bash
	$ INPUT_JSON=https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK_subsampled_chr19_only_caper.json
	```

2-1) **Conda**: Run a workflow with Conda. Make sure that you have followed [this instruction](docs/install_conda.md) to install Conda and its environments.
	> **WARNING**: We no longer recommend Conda for resolving dependencies and plan to phase out Conda support. Instead we recommend using Docker or Singularity. You can install Singularity and use it for our pipeline with Caper (by adding `--singularity` to command line arguments).

	```bash
	$ conda activate encode-chip-seq-pipeline
	$ caper run chip.wdl -i ${INPUT_JSON}
	```

2-2) **Singularity (RECOMMENDED)**: Run a workflow with Singularity.

	```bash
	$ caper run chip.wdl -i ${INPUT_JSON} --singularity
	```

	> **HPCs**: To run multiple workflows on HPCs (e.g. Stanford Sherlock and SCG) see details at [Caper's README](https://github.com/ENCODE-DCC/caper/blob/master/README.md#how-to-run-it-on-slurm-cluster). Do not run Caper on login nodes. Your workflows will get killed. There is a learning curve to understand server/client structure of Caper/Cromwell.


2-3) **Docker**: Run a workflow with Docker.

  ```bash
	$ caper run chip.wdl -i ${INPUT_JSON} --docker
  ```

3) You can also run a workflow on cloud platforms such as AWS (`aws`) and Google Cloud Platform (`gcp`) if Caper's configuration file is correctly configured for them. See details at [Caper's README](https://github.com/ENCODE-DCC/caper).


## Running pipelines without Caper

Caper uses the cromwell workflow execution engine to run the workflow on the platform you specify.  While we recommend you use caper, if you want to run cromwell directly without caper you can learn about that [here](docs/deprecated/OLD_METHOD.md).

## DNAnexus

You can also run our pipeline on DNAnexus without using Caper or Cromwell. There are two ways to build a workflow on DNAnexus based on our WDL.

1) [dxWDL CLI](docs/tutorial_dx_cli.md)
2) [DNAnexus Web UI](docs/tutorial_dx_web.md)

## Input JSON file

An input JSON file includes all genomic data files, input parameters and metadata for running pipelines. Always use absolute paths in an input JSON.

[Input JSON file specification](docs/input.md)

## How to organize outputs

Install [Croo](https://github.com/ENCODE-DCC/croo#installation). Make sure that you have python3(> 3.4.1) installed on your system. Find a `metadata.json` on Caper's output directory.

```bash
$ pip install croo
$ croo [METADATA_JSON_FILE]
```

## Useful tools

There are some useful tools to post-process outputs of the pipeline.

### qc_jsons_to_tsv

[This tool](utils/qc_jsons_to_tsv/README.md) recursively finds and parses all `qc.json` (pipeline's [final output](docs/example_output/v1.1.5/qc.json)) found from a specified root directory. It generates a TSV file that has all quality metrics tabulated in rows for each experiment and replicate. This tool also estimates overall quality of a sample by [a criteria definition JSON file](utils/qc_jsons_to_tsv/criteria.default.json) which can be a good guideline for QC'ing experiments.
