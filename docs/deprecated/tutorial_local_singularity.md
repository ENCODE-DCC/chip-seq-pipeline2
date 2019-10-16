# Tutorial for general UNIX computers with singularity

1. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ cd
    $ wget https://github.com/broadinstitute/cromwell/releases/download/38/cromwell-38.jar
    $ chmod +rx cromwell-38.jar
    ```

2. Git clone this pipeline and move into it.
    ```bash
    $ cd
    $ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
    $ cd chip-seq-pipeline2
    ```

3. Download a SUBSAMPLED paired-end sample of [ENCSR936XTK](https://www.encodeproject.org/experiments/ENCSR936XTK/).
    ```bash
    $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/ENCSR936XTK_fastq_subsampled.tar
    $ tar xvf ENCSR936XTK_fastq_subsampled.tar
    ```

4. Download pre-built chr19/chrM-only genome database for hg38.
    ```bash
    $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_chr19_chrM_chip.tar
    $ tar xvf test_genome_database_hg38_chr19_chrM_chip.tar
    ```

5. CHECK YOUR SINGULARITY VERSION FIRST AND UPGRADE IT TO A VERSION `>=2.5.2` OR PIPELINE WILL NOT WORK CORRECTLY.
    ```bash
    $ singularity --version
    ```

6. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```bash
    $ mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name chip-seq-pipeline-dev-v1.3.2.simg -F docker://quay.io/encode-dcc/chip-seq-pipeline:dev-v1.3.2
    ```

7. Run a pipeline for the test sample.
    ```bash
    $ INPUT=dev/examples/local/ENCSR936XTK_subsampled_chr19_only.json
    $ PIPELINE_METADATA=metadata.json
    $ java -jar -Xmx1G -Dconfig.file=dev/backends/backend.conf -Dbackend.default=singularity cromwell-38.jar run chip.wdl -i ${INPUT} -o dev/workflow_opts/singularity.json -m ${PIPELINE_METADATA}
    ```

8. It will take about 6 hours. You will be able to find all outputs on `cromwell-executions/chip/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

10. You can resume a failed pipeline from where it left off by using `PIPELINE_METADATA`(`metadata.json`) file. This file is created for each pipeline run. See [here](../utils/resumer/README.md) for details. Once you get a new input JSON file from the resumer, use it `INPUT=resume.[FAILED_WORKFLOW_ID].json` instead of `INPUT=dev/examples/local/ENCSR936XTK_subsampled_chr19_only.json`.

11. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `dev/workflow_opts/singularity.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/chip-seq-pipeline-dev-v1.3.2.simg",
            "singularity_bindpath" : "/your/,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR2,..."
        }
    }
    ```
