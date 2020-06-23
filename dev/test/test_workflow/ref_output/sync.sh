#!/bin/bash

gsutil -m rsync -r -d . gs://encode-pipeline-test-samples/encode-chip-seq-pipeline/ref_output
