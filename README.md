# ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline

[![CircleCI](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master)


## Important notice for Conda users

If it takes too long to resolve Conda package conflicts while installing pipeline's Conda environment, then try with `mamba` instead. Add `mamba` to the install command line.
```bash
$ scripts/install_conda_env.sh mamba
```

For every new pipeline release, Conda users always need to update pipeline's Conda environment (`encode-chip-seq-pipeline`) even though they don't use new added features.
```bash
$ cd chip-seq-pipeline2
$ scripts/update_conda_env.sh
```

For pipelines >= v1.7.0, Conda users also need to manually install Caper and Croo **INSIDE** the environment. These two tools have been removed from pipeline's Conda environment since v1.7.0.
```bash
$ source activate encode-chip-seq-pipeline
$ pip install caper croo
```


## Redacting filtered/deduped BAM (new feature >= v1.7.0)

Added `"chip.redact_nodup_bam"` to redact (de-identify) BAM. Such conversion is done with `ptools` which is not officially registered to PIP and Conda repository (`bioconda`) yet.

> **IMPORTANT**: Alignment quality metrics calculated during/before filtering (task `filter`) will still be based on non-redacted original BAMs. However, all downstream analyses (e.g. peak-calling) will be based on redact nodup BAM.

Conda users need to install it from our temporary PIP repo (`ptools_bin`). GCP/AWS/Docker/Singularity users do not need to install it. `ptools` is already included in new pipeline's Docker/Singularity image. This is for Conda users only. We will remove this warning in the next release after `ptools` is registered to `bioconda` and added to the requirements list (`scripts/requirements.txt`).

```
# Activate pipeline's conda environment
$ source activate encode-chip-seq-pipeline

# Install ptools_bin inside the environment
(encode-chip-seq-pipeline) $ pip3 install ptools_bin

# Run pipelines in the environment
(encode-chip-seq-pipeline) $ caper run/submit ...
```

## Introduction 
This ChIP-Seq pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications (by Anshul Kundaje) in [this google doc](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#).

### Features

* **Portability**: The pipeline run can be performed across different cloud platforms such as Google, AWS and DNAnexus, as well as on cluster engines such as SLURM, SGE and PBS.
* **User-friendly HTML report**: In addition to the standard outputs, the pipeline generates an HTML report that consists of a tabular representation of quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures). An example of the [HTML report](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/example_output/qc.html). The [json file](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/example_output/qc.json) used in generating this report.
* **Supported genomes**: Pipeline needs genome specific data such as aligner indices, chromosome sizes file and blacklist. We provide a genome database downloader/builder for hg38, hg19, mm10, mm9. You can also use this [builder](docs/build_genome_database.md) to build genome database from FASTA for your custom genome.

## Installation

1) Git clone this pipeline.
	> **IMPORTANT**: use `~/chip-seq-pipeline2/chip.wdl` as `[WDL]` in Caper's documentation.

	```bash
	$ cd
	$ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
	```

2) Install pipeline's [Conda environment](docs/install_conda.md) if you want to use Conda instead of Docker/Singularity. Conda is recommneded on local computer and HPCs (e.g. Stanford Sherlock/SCG). Use 
	> **IMPORTANT**: use `encode-chip-seq-pipeline` as `[PIPELINE_CONDA_ENV]` in Caper's documentation.
 
3) **Skip this step if you have installed pipeline's Conda environment**. Caper is already included in the Conda environment. [Install Caper](https://github.com/ENCODE-DCC/caper#installation). Caper is a python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell).

	> **IMPORTANT**: Make sure that you have python3(>= 3.6.0) installed on your system.

	```bash
	$ pip install caper  # use pip3 if it doesn't work
	```

4) Follow [Caper's README](https://github.com/ENCODE-DCC/caper) carefully. Find an instruction for your platform.
	> **IMPORTANT**: Configure your Caper configuration file `~/.caper/default.conf` correctly for your platform.

## Test input JSON file

Use `https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI_subsampled_chr19_only.json` as `[INPUT_JSON]` in Caper's documentation.

## Input JSON file

> **IMPORTANT**: DO NOT BLINDLY USE A TEMPLATE/EXAMPLE INPUT JSON. READ THROUGH THE FOLLOWING GUIDE TO MAKE A CORRECT INPUT JSON FILE.

An input JSON file specifies all the input parameters and files that are necessary for successfully running this pipeline. This includes a specification of the path to the genome reference files and the raw data fastq file. Please make sure to specify absolute paths rather than relative paths in your input JSON files.

1) [Input JSON file specification (short)](docs/input_short.md)
2) [Input JSON file specification (long)](docs/input.md)

## Running and sharing on Truwl
You can run this pipeline on [truwl.com](https://truwl.com/). This provides a web interface that allows you to define inputs and parameters, run the job on GCP, and monitor progress. To run it you will need to create an account on the platform then request early access by emailing [info@truwl.com](mailto:info@truwl.com) to get the right permissions. You can see the example cases from this repo at [https://truwl.com/workflows/instance/WF_dd6938.8f.340f/command](https://truwl.com/workflows/instance/WF_dd6938.8f.340f/command) and [https://truwl.com/workflows/instance/WF_dd6938.8f.8aa3/command](https://truwl.com/workflows/instance/WF_dd6938.8f.8aa3/command). The example jobs (or other jobs) can be forked to pre-populate the inputs for your own job.

If you do not run the pipeline on Truwl, you can still share your use-case/job on the platform by getting in touch at [info@truwl.com](mailto:info@truwl.com) and providing your inputs.json file.

## Running a pipeline on DNAnexus

You can also run this pipeline on DNAnexus without using Caper or Cromwell. There are two ways to build a workflow on DNAnexus based on our WDL.

1) [dxWDL CLI](docs/tutorial_dx_cli.md)
2) [DNAnexus Web UI](docs/tutorial_dx_web.md)

## How to organize outputs

Install [Croo](https://github.com/ENCODE-DCC/croo#installation). **You can skip this installation if you have installed pipeline's Conda environment and activated it**. Make sure that you have python3(> 3.4.1) installed on your system. Find a `metadata.json` on Caper's output directory.

```bash
$ pip install croo
$ croo [METADATA_JSON_FILE]
```

## How to make a spreadsheet of QC metrics

Install [qc2tsv](https://github.com/ENCODE-DCC/qc2tsv#installation). Make sure that you have python3(> 3.4.1) installed on your system. 

Once you have [organized output with Croo](#how-to-organize-outputs), you will be able to find pipeline's final output file `qc/qc.json` which has all QC metrics in it. Simply feed `qc2tsv` with multiple `qc.json` files. It can take various URIs like local path, `gs://` and `s3://`.

```bash
$ pip install qc2tsv
$ qc2tsv /sample1/qc.json gs://sample2/qc.json s3://sample3/qc.json ... > spreadsheet.tsv
```

QC metrics for each experiment (`qc.json`) will be split into multiple rows (1 for overall experiment + 1 for each bio replicate) in a spreadsheet.


## Troubleshooting

See [this document](docs/troubleshooting.md) for troubleshooting. I will keep updating this document for errors reported by users.
