Dev
===

## Building templates on DX for each genome
```
# general
java -jar ~/dxWDL-0.75.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/v1.1/general -defaults examples/dx/template_general.json

# hg38
java -jar ~/dxWDL-0.75.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/v1.1/hg38 -defaults examples/dx/template_hg38.json

# hg19
java -jar ~/dxWDL-0.75.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/v1.1/hg19 -defaults examples/dx/template_hg19.json

# mm10
java -jar ~/dxWDL-0.75.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/v1.1/mm10 -defaults examples/dx/template_mm10.json

# mm9
java -jar ~/dxWDL-0.75.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/v1.1/mm9 -defaults examples/dx/template_mm9.json

# test sample
java -jar ~/dxWDL-0.75.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/v1.1/test_ENCSR936XTK_subsampled -defaults examples/dx/ENCSR936XTK_subsampled_dx.json
```
