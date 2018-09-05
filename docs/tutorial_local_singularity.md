Tutorial for general UNIX computers with singularity
====================================================

1. Git clone this pipeline and move into it.
    ```
      $ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
      $ cd chip-seq-pipeline2
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```
      $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
      $ chmod +rx cromwell-34.jar
    ```

3. Download a SUBSAMPLED paired-end sample of [ENCSR936XTK](https://www.encodeproject.org/experiments/ENCSR936XTK/).
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/ENCSR936XTK_fastq_subsampled.tar
      $ tar xvf ENCSR936XTK_fastq_subsampled.tar
    ```

4. Download pre-built genome database for hg38.
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_chip.tar
      $ tar xvf test_genome_database_hg38_chip.tar
    ```

5. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```
      $ SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://quay.io/encode-dcc/chip-seq-pipeline:v1.1
    ```

6. [Bind your input/genome data directories to singularity](singularity.md). You can skip this step for this example. However, if you want to use your own input data and genome database then you may need to read through it very carefully and ADD YOUR DIRECTORIES TO `workflow_opts/singularity.json`.

7. Run a pipeline for the test sample.
    ```
      $ INPUT=examples/local/ENCSR936XTK_subsampled.json
      $ java -jar -Dconfig.file=backends/backend.conf cromwell-34.jar run chip.wdl -i ${INPUT} -o workflow_opts/singularity.json
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/chip/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

10. Please read through the next section. It's VERY IMPORTANT.
