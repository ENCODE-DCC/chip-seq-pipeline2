WORKDIR=/srv/scratch/shared/surya/leepc12/run/chip-seq-pipeline/cromwell/ENCSR000DYI_chr21/from_fastq
mkdir -p $WORKDIR && cd $WORKDIR
source activate encode-chip-seq-pipeline
java -jar -Dconfig.file=$CODE/chip-seq-pipeline-dev/backends/default.conf $CODE/chip-seq-pipeline-dev/cromwell-30-x.jar run $CODE/chip-seq-pipeline-dev/chipseq.wdl -i $CODE/chip-seq-pipeline-dev/examples/klab/ENCSR000DYI_chr21_klab.json -o $CODE/chip-seq-pipeline-dev/workflow_opts/non_docker.json -m output_ENCSR000DYI.json
source deactivate
