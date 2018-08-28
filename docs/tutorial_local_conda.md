Tutorial for general UNIX computers without docker
==================================================

1. Git clone this pipeline.
    ```
      $ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
    ```

2. Move to pipeline's directory.
    ```
      $ cd chip-seq-pipeline2
    ```

3. [Install Conda](https://conda.io/miniconda.html)

4. Install Conda dependencies.
    ```
      $ bash installers/uninstall_dependencies.sh  # to remove any existing pipeline env
      $ bash installers/install_dependencies.sh
    ```

5. Download cromwell.
    ```
      $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
      $ chmod +rx cromwell-34.jar
    ```

6. Download a SUBSAMPLED paired-end sample of [ENCSR936XTK](https://www.encodeproject.org/experiments/ENCSR936XTK/).
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/ENCSR936XTK_fastq_subsampled.tar
      $ tar xvf ENCSR936XTK_fastq_subsampled.tar
    ```

7. Download pre-built genome database for hg38.
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_chip.tar
      $ tar xvf test_genome_database_hg38_chip.tar
    ```

8. Run a pipeline for the test sample.
    ```
      $ source activate encode-chip-seq-pipeline # IMPORTANT!
      $ INPUT=examples/local/ENCSR936XTK_subsampled.json
      $ java -jar -Dconfig.file=backends/backend.conf cromwell-34.jar run chip.wdl -i ${INPUT}
    ```

9. It will take about an hour. You will be able to find all outputs on `cromwell-executions/chip/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

10. See full specification for [input JSON file](input.md).
