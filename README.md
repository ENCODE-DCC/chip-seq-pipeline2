# ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline

[![CircleCI](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/chip-seq-pipeline2/tree/master)


## Introduction

This ChIP-Seq pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications (by Anshul Kundaje) in [this google doc](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#).

### Features

* **Portability**: The pipeline run can be performed across different cloud platforms such as Google, AWS and DNAnexus, as well as on cluster engines such as SLURM, SGE and PBS.
* **User-friendly HTML report**: In addition to the standard outputs, the pipeline generates an HTML report that consists of a tabular representation of quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures). An example of the [HTML report](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/example_output/qc.html). The [json file](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/example_output/qc.json) used in generating this report.
* **Supported genomes**: Pipeline needs genome specific data such as aligner indices, chromosome sizes file and blacklist. We provide a genome database downloader/builder for hg38, hg19, mm10, mm9. You can also use this [builder](docs/build_genome_database.md) to build genome database from FASTA for your custom genome.

## Installation

1) Install Caper (Python Wrapper/CLI for [Cromwell](https://github.com/broadinstitute/cromwell)).
	```bash
	$ pip install caper
	```

2) **IMPORTANT**: Read Caper's [README](https://github.com/ENCODE-DCC/caper/blob/master/README.md) carefully to choose a backend for your system. Follow the instruction in the configuration file.
	```bash
	# backend: local or your HPC type (e.g. slurm, sge, pbs, lsf). read Caper's README carefully.
	$ caper init [YOUR_BACKEND]

	# IMPORTANT: edit the conf file and follow commented instructions in there
	$ vi ~/.caper/default.conf
	```

3) Git clone this pipeline.
	```bash
	$ cd
	$ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
	```

4) Define test input JSON.
	```bash
	INPUT_JSON="https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI_subsampled_chr19_only.json"
	```

5) If you have Docker and want to run pipelines locally on your laptop. `--max-concurrent-tasks 1` is to limit number of concurrent tasks to test-run the pipeline on a laptop. Uncomment it if run it on a workstation/HPC.
	```bash
	# check if Docker works on your machine
	$ docker run ubuntu:latest echo hello

	# --max-concurrent-tasks 1 is for computers with limited resources
	$ caper run chip.wdl -i "${INPUT_JSON}" --docker --max-concurrent-tasks 1
	```

6) Otherwise, install Singularity on your system. Please follow [this instruction](https://neuro.debian.net/install_pkg.html?p=singularity-container) to install Singularity on a Debian-based OS. Or ask your system administrator to install Singularity on your HPC.
	```bash
	# check if Singularity works on your machine
	$ singularity exec docker://ubuntu:latest echo hello

	# on your local machine (--max-concurrent-tasks 1 is for computers with limited resources)
	$ caper run chip.wdl -i "${INPUT_JSON}" --singularity --max-concurrent-tasks 1

	# on HPC, make sure that Caper's conf ~/.caper/default.conf is correctly configured to work with your HPC
    # the following command will submit Caper as a leader job to SLURM with Singularity
    $ caper hpc submit chip.wdl -i "${INPUT_JSON}" --singularity --leader-job-name ANY_GOOD_LEADER_JOB_NAME

    # check job ID and status of your leader jobs
    $ caper hpc list

    # cancel the leader node to close all of its children jobs
    # If you directly use cluster command like scancel or qdel then
    # child jobs will not be terminated
    $ caper hpc abort [JOB_ID]
	```

7) (Optional Conda method) **WE DO NOT HELP USERS FIX CONDA DEPENDENCY ISSUES. IF CONDA METHOD FAILS THEN PLEASE USE SINGULARITY METHOD INSTEAD**. **DO NOT USE A SHARED CONDA. INSTALL YOUR OWN [MINICONDA3](https://docs.conda.io/en/latest/miniconda.html) AND USE IT.**
	```bash
	# check if you are not using a shared conda, if so then delete it or remove it from your PATH
	$ which conda

	# uninstall pipeline's old environments
	$ bash scripts/uninstall_conda_env.sh

	# install new envs, you need to run this for every pipeline version update.
	# it may be killed if you run this command line on a login node on HPC.
	# it's recommended to make an interactive node with enough resources and run it there.
	$ bash scripts/install_conda_env.sh

	# if installation fails please use Singularity method instead.

	# on your local machine (--max-concurrent-tasks 1 is for computers with limited resources)
	$ caper run chip.wdl -i "${INPUT_JSON}" --conda --max-concurrent-tasks 1

	# on HPC, make sure that Caper's conf ~/.caper/default.conf is correctly configured to work with your HPC
    # the following command will submit Caper as a leader job to SLURM with Conda
    $ caper hpc submit chip.wdl -i "${INPUT_JSON}" --conda --leader-job-name ANY_GOOD_LEADER_JOB_NAME

    # check job ID and status of your leader jobs
    $ caper hpc list

    # cancel the leader node to close all of its children jobs
    # If you directly use cluster command like scancel or qdel then
    # child jobs will not be terminated
    $ caper hpc abort [JOB_ID]
	```


## Input JSON file

> **IMPORTANT**: DO NOT BLINDLY USE A TEMPLATE/EXAMPLE INPUT JSON. READ THROUGH THE FOLLOWING GUIDE TO MAKE A CORRECT INPUT JSON FILE.

An input JSON file specifies all the input parameters and files that are necessary for successfully running this pipeline. This includes a specification of the path to the genome reference files and the raw data fastq file. Please make sure to specify absolute paths rather than relative paths in your input JSON files.

1) [Input JSON file specification (short)](docs/input_short.md)
2) [Input JSON file specification (long)](docs/input.md)


## Running on Terra/Anvil (using Dockstore)

Visit our pipeline repo on [Dockstore](https://dockstore.org/workflows/github.com/ENCODE-DCC/chip-seq-pipeline2). Click on `Terra` or `Anvil`. Follow Terra's instruction to create a workspace on Terra and add Terra's billing bot to your Google Cloud account.

Download this [test input JSON for Terra](https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI_subsampled_chr19_only.terra.json) and upload it to Terra's UI and then run analysis.

If you want to use your own input JSON file, then make sure that all files in the input JSON are on a Google Cloud Storage bucket (`gs://`). URLs will not work.


## Running on DNAnexus (using Dockstore)

Sign up for a new account on [DNAnexus](https://platform.dnanexus.com/) and create a new project on either AWS or Azure. Visit our pipeline repo on [Dockstore](https://dockstore.org/workflows/github.com/ENCODE-DCC/chip-seq-pipeline2). Click on `DNAnexus`. Choose a destination directory on your DNAnexus project. Click on `Submit` and visit DNAnexus. This will submit a conversion job so that you can check status of it on `Monitor` on DNAnexus UI.

Once conversion is done download one of the following input JSON files according to your chosen platform (AWS or Azure) for your DNAnexus project:
- AWS: https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI_subsampled_chr19_only_dx.json
- Azure: https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI_subsampled_chr19_only_dx_azure.json

You cannot use these input JSON files directly. Go to the destination directory on DNAnexus and click on the converted workflow `chip`. You will see input file boxes in the left-hand side of the task graph. Expand it and define FASTQs (`fastq_repX_R1` and `fastq_repX_R1`) and `genome_tsv` as in the downloaded input JSON file. Click on the `common` task box and define other non-file pipeline parameters.  e.g. `pipeline_type`, `paired_end` and `ctl_paired_end`.

We have a separate project on DNANexus to provide example FASTQs and `genome_tsv` for `hg38` and `mm10` (also chr19-only version of those two. Use chr19-only versions for testing). We recommend to make copies of these directories on your own project.

`genome_tsv`
- AWS: https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v4
- Azure: https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v4

Example FASTQs
- AWS: https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/fastq_subsampled
- Azure: https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-test-samples/encode-chip-seq-pipeline/ENCSR000DYI/fastq_subsampled


## Running on DNAnexus (using our pre-built workflows)

See [this](docs/tutorial_dx_web.md) for details.


## Running and sharing on Truwl
You can run this pipeline on [truwl.com](https://truwl.com/). This provides a web interface that allows you to define inputs and parameters, run the job on GCP, and monitor progress. To run it you will need to create an account on the platform then request early access by emailing [info@truwl.com](mailto:info@truwl.com) to get the right permissions. You can see the example cases from this repo at [https://truwl.com/workflows/instance/WF_dd6938.8f.340f/command](https://truwl.com/workflows/instance/WF_dd6938.8f.340f/command) and [https://truwl.com/workflows/instance/WF_dd6938.8f.8aa3/command](https://truwl.com/workflows/instance/WF_dd6938.8f.8aa3/command). The example jobs (or other jobs) can be forked to pre-populate the inputs for your own job.

If you do not run the pipeline on Truwl, you can still share your use-case/job on the platform by getting in touch at [info@truwl.com](mailto:info@truwl.com) and providing your inputs.json file.

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
