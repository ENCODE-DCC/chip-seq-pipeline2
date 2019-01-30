# Dev

## Command line for version change
```bash
PREV_VER=v1.1.5
NEW_VER=v1.1.6
for f in $(grep -rl ${PREV_VER} --include=*.{wdl,md,sh,yml})
do
  sed -i "s/${PREV_VER}/${NEW_VER}/g" ${f}
done
cd workflow_opts
for f in $(grep -rl ${PREV_VER} --include=*.json)
do
  sed -i "s/${PREV_VER}/${NEW_VER}/g" ${f}
done
cd ..
```

## Building templates on DX for each genome

Make sure that you have [`dxWDL-0.77.jar`](https://github.com/DNAnexus/dxWDL/releases/download/0.77/dxWDL-0.77.jar) on your `$HOME`. Install [DNAnexus Platform SDK](https://wiki.DNAnexus.com/downloads) with `pip install dxpy`. Log-in on DNAnexus with `dx login` and choose "ENCODE Uniform Processing Pipelines" (name of our official DNAnexus project for pipelines).

Run the following command line locally to build out DX workflows for this pipeline on our official one. This will overwrite (`-f` parameter does it).

```bash
# version
VER=v1.1.6

# general
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/general -defaults examples/dx/template_general.json

# hg38
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/hg38 -defaults examples/dx/template_hg38.json

# hg19
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/hg19 -defaults examples/dx/template_hg19.json

# mm10
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/mm10 -defaults examples/dx/template_mm10.json

# mm9
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/mm9 -defaults examples/dx/template_mm9.json

# test sample PE ENCSR936XTK (full)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR936XTK -defaults examples/dx/ENCSR936XTK_dx.json

# test sample SE ENCSR000DYI (full)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR000DYI -defaults examples/dx/ENCSR000DYI_dx.json

# test sample PE ENCSR936XTK (subsampled, chr19/chrM only)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR936XTK_subsampled_chr19_only -defaults examples/dx/ENCSR936XTK_subsampled_chr19_only_dx.json

# test sample SE ENCSR000DYI (subsampled, chr19/chrM only)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR000DYI_subsampled_chr19_only -defaults examples/dx/ENCSR000DYI_subsampled_chr19_only_dx.json

## DX Azure

# general
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/general -defaults examples/dx_azure/template_general.json

# hg38
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/hg38 -defaults examples/dx_azure/template_hg38.json

# hg19
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/hg19 -defaults examples/dx_azure/template_hg19.json

# mm10
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/mm10 -defaults examples/dx_azure/template_mm10.json

# mm9
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/mm9 -defaults examples/dx_azure/template_mm9.json

# test sample PE ENCSR936XTK (full)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR936XTK -defaults examples/dx_azure/ENCSR936XTK_dx_azure.json

# test sample SE ENCSR000DYI (full)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR000DYI -defaults examples/dx_azure/ENCSR000DYI_dx_azure.json

# test sample PE ENCSR936XTK (subsampled, chr19/chrM only)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR936XTK_subsampled_chr19_only -defaults examples/dx_azure/ENCSR936XTK_subsampled_chr19_only_dx_azure.json

# test sample SE ENCSR000DYI (subsampled, chr19/chrM only)
java -jar ~/dxWDL-0.77.jar compile chip.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ChIP-seq2/workflows/$VER/test_ENCSR000DYI_subsampled_chr19_only -defaults examples/dx_azure/ENCSR000DYI_subsampled_chr19_only_dx_azure.json
```
