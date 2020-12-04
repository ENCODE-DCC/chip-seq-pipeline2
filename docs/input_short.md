# Input JSON

An input JSON file is a file which must include all the information needed to run this pipeline. Hence, it must include the absolute paths to all the control and experimental fastq files; paths to all the genomic data files needed for this pipeline, and it must also specify the parameters and the metadata needed for running this pipeline. If the parameters are not specified in an input JSON file, default values will be used. We provide a set of template JSON files: [minimum](../example_input_json/template.json) and [full](../example_input_json/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

>**IMPORTANT**: ALWAYS USE ABSOLUTE PATHS.

# Checklist

Mandatory parameters.

1) Pipeline type
    * `chip.pipeline_type`: `tf` for TF ChIP-seq, `histone` for histone ChIP-seq and `control` for aligning control FASTQs only. One major difference between two types is that `tf` uses `spp` peak caller with controls but `histone` uses `macs2` peak caller without controls. If it is `control` then `chip.align_only` is automatically turned on and pipeline will align FASTQs defined in `chip.fastqs_repX_RY` (where `X` means control replicate ID and `Y` is read-end) as controls. For control mode, do not define  `chip.ctl_fastqs_repX_RY`. Pipeline will skip cross-correlation/JSD/GC-bias analysis for control mode.

2) Experiment title/description
    * `chip.title`: experiment title for a final HTML report.
    * `chip.description`: experiment description for a final HTML report.

3) Read endedness
    * `chip.paired_end`: `true` if ALL replicates are paired-ended.
    * (Optional) `chip.paired_ends`: For samples with mixed read ends, you can define read endedness for each biological replicate (e.g. `[true, false]` means paired-ended biorep-1 and single-ended biorep-2).
    * `chip.ctl_paired_end`: `true` if ALL controls are paired-ended. If not defined then `chip.paired_end` will be used.
    * (Optional) `chip.ctl_paired_ends`: For controls with mixed read ends, you can define read endedness for each biological replicate (e.g. `[true, false]` means paired-ended biorep-1 and single-ended biorep-2). If not defined then `chip.paired_ends` will be used.

4) Reference genome
    * `chip.genome_tsv`: Choose one from the following genome TSVs. `v3` is a standard for >=ENCODE4.
        Genome|URL
        -|-
        hg38|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv`
        mm10|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/mm10.tsv`
        hg19|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/hg19_caper.tsv`
        mm9|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/mm9_caper.tsv`

        For DNAnexus CLI (AWS project):
        Genome|DX URI
        -|-
        hg38|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v3/hg38.dx.tsv`
        mm10|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v3/mm10.dx.tsv`

        For DNAnexus CLI (Azure project): 
        Genome|DX URI
        -|-
        hg38|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v3/hg38.dx_azure.tsv`
        mm10|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v3/mm10.dx_azure.tsv`

        For DNAnexus Web UI (AWS project): Choose one of the following TSV file on `https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v3`.
        Genome|File name
        -|-
        hg38|`hg38.dx.tsv`
        mm10|`mm10.dx.tsv`

        For DNAnexus Web UI (Azure project): Choose one of the following TSV file on `https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v3`.
        Genome|File name
        -|-
        hg38|`hg38.dx_azure.tsv`
        mm10|`mm10.dx_azure.tsv`

    * To build a new TSV file from use your own FASTA (`.fa` and `.2bit`) see [this](build_genome_database.md).
    * You can also define genome specific parameters (defined in a genome TSV file) in an input JSON file. Parameters defined in an input JSON file will override those defined in a genome TSV file. For example, you can simply replace a blacklist file while keeping all other parameters in a genome TSV file for `hg38`.
        ```javascript
        {
            "chip.genome_tsv" : "/somewhere/hg38.tsv",
            "chip.blacklist": "/new/genome/data/new_blacklist.bed.gz"
        }
        ```

5) [Input files](#input-files)
    * See [this](#input-files) for how to define FASTQ/BAM/TAG-ALIGNs for your sample.
    * ChIP-seq pipeline does not have an adapter detector/trimmer.

6) Important parameters
    * `chip.crop_length`: Crop FASTQs with Trimmomatic (using parameters `CROP`). It is 0 (disabled) by default.
    * `chip.crop_length_tol`: Read length tolerance while cropping FASTQs with `chip.crop_length`. Trimmomatic's `MINLEN` will be set as `chip.crop_length` - `abs(chip.crop_length_tol)`. **WARNING**: Check your FASTQs' read length first. Reads SHORTER than `MINLEN` will be excluded while cropping, hence not included in output BAM files and all downstream analyses. Tolerance is 2 by default.
    * `chip.pval_thresh`: P-value threshold for peak-caller MACS2 (macs2 callpeak -p).
    * `chip.fdr_thresh`: FDR threshold for peak-caller SPP (run_spp.R -fdr=).
    * `chip.idr_thresh`: IDR (irreproducible discovery rate) threshold.

7) Parameters for subsampling experiment replicates and controls.
    These parameters are used for subsampling on filtered (nodup) BAMs. Therefore, setting these parameters will affect all downsteam analyses including peak-calling. It's 0 by default, which means no subsampling. These parameters will be applied to each experiment replicate and control.
   
    * `chip.subsample_reads`: Subsample experiment replicate's reads. For PE dataset, this is not a number of read pairs but number of reads. 
    * `chip.ctl_subsample_reads`: Subsample control's reads. For PE dataset, this is not a number of read pairs but number of reads. 

8) Parameters for automatically choosing an appropriate control and subsampling it.
    Before calling peaks, there is a separate (from item 7) mechanism to subsample controls to prevent using too deep controls. Pipeline has a logic to find an approriate control for each experiment replicate. Choices are 1) control with the same index (e.g. rep2 vs. ctl2), 2) pooled control (e.g. rep2 vs. ctl1 + ctl2).
    * `chip.always_use_pooled_ctl`: (For experiments with controls only) Always use a pooled control to compare with each replicate. If a single control is given then use it. It is enabled by default.
    * `chip.ctl_depth_ratio`: (For experiments with controls only) If ratio of depth between controls is higher than this. then always use a pooled control for all replicates. It's 1.2 by default.

    > **IMPORTANT**: Either of the following parameters are set to > 0 for automatic subsampling. Set all of them to 0 to disable automatic subsample. Such automatic control subsampling does not work with controls with mixed endedness (e.g. SE ctl-rep1 + PE ctl-rep2). Pipeline calculates `Limit1` and `Limit2` from the following two parameters. And then takes a maximum of them and calls it `Limit`. If a chosen control for experiment replicate (each one and pooled one) is deeper than `Limit` then that such control is subsampled to `Limit`.
    * `chip.exp_ctl_depth_ratio_limit`: (For experiments with controls only) `Limit1` is defined as experiment replicate's depth multiplied by this parameter. It's 5.0 by default.
    * `chip.ctl_depth_limit`: (For experiments with controls only) Hard limit on control's depth. If control is deeper than this hard limit then such control is subsampled to it.

9) [Resources](#resources)
    * It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task.

Optional parameters.

10) Other useful parameters
    * `chip.fraglen`: Array of Integers. Fragment length for each bio replicate. If you start from FASTQs then our pipeline automatically estimate it from cross-correlation analysis (task `xcor`) result since such analysis requires a special treamtment for FASTQs. It is possible that fragment length is not estimated correctly (or pipeline can fail due to negative fraglen) if you start from different types (BAM/TAG-ALIGN). For such case, you can manually define fragment length for each bio rep. (e.g. `[200, 150]` means 200 for rep1 and 150 for rep2).

11) Flags
    * `chip.redact_bam`: BAMs aligned from FASTQs will be redacted. Variations like indels will be replaced with corresponding reference sequence. **THIS FLAG IS FOR PIPELINE STARTING FROM FASTQS ONLY**. Starting from BAMs will deactivate this flag.
    * `chip.align_only`: Peak calling and its downstream analyses will be disabled. Useful if you just want to align your FASTQs into filtered BAMs/TAG-ALIGNs and don't want to call peaks on them. **IMPORTANT**: Even though `chip.pipeline_type` does not matter for align only mode, you still need to define it since it is a required parameter in WDL. Define it as `tf` for such cases.
    * `chip.true_rep_only`: Disable pseudo replicate generation and all related analyses


## Input files

> **IMPORTANT**: Our pipeline considers a replicate (`rep`) as a biological replicate. You can still define technical replicates for each bio replicate. Tech replicates will be merged together to make a single FASTQ for each bio replicate. Controls can also have technical replicates.

> **IMPORTANT**: Our pipeline supports up to 10 bio replicates and 10 controls.

> **IMPORTANT**: Our pipeline has cross-validation analyses (IDR/overlap) comparing every pair of all replicates. Number of tasks for such analyses will be like <sub>n</sub>C<sub>2</sub>. This number will be 45 for 10 bio replicates. It's recommended to keep number of replicates <= 4.

Pipeline can start from any of the following data types (FASTQ, BAM, NODUP_BAM and TAG-ALIGN). 

1) Starting from FASTQs
    * Technical replicates for each bio-rep will be **MERGED** in the very early stage of the pipeline. Each read end R1 and R2 have separate arrays `chip.fastqs_repX_R1` and `chip.fastqs_repX_R2`. Do not define R2 array for single-ended replicates.
    * Example of 3 paired-ended biological replicates and 2 technical replicates for each bio rep. Two technical replicates `BIOREPX_TECHREP1.R1.fq.gz` and `BIOREPX_TECHREP2.R1.fq.gz` for each bio replicate will be merged.

        ```javascript
        {
            "chip.paired_end" : true,
            "chip.fastqs_rep1_R1" : ["BIOREP1_TECHREP1.R1.fq.gz", "BIOREP1_TECHREP2.R1.fq.gz"],
            "chip.fastqs_rep1_R2" : ["BIOREP1_TECHREP1.R2.fq.gz", "BIOREP1_TECHREP2.R2.fq.gz"],
            "chip.fastqs_rep2_R1" : ["BIOREP2_TECHREP1.R1.fq.gz", "BIOREP2_TECHREP2.R1.fq.gz"],
            "chip.fastqs_rep2_R2" : ["BIOREP2_TECHREP1.R2.fq.gz", "BIOREP2_TECHREP2.R2.fq.gz"],
            "chip.fastqs_rep3_R1" : ["BIOREP3_TECHREP1.R1.fq.gz", "BIOREP3_TECHREP2.R1.fq.gz"],
            "chip.fastqs_rep3_R2" : ["BIOREP3_TECHREP1.R2.fq.gz", "BIOREP3_TECHREP2.R2.fq.gz"]
        }
        ```

2) Starting from BAMs
    * Define a BAM for each replicate. Our pipeline does not determine read endedness from a BAM file. You need to explicitly define read endedness.
    * Example of 3 singled-ended replicates.
        ```javascript
        {       
            "chip.paired_end" : false,
            "chip.bams" : ["rep1.bam", "rep2.bam", "rep3.bam"]
        }
        ```

3) Starting from filtered/deduped BAMs
    * Define a filtered/deduped BAM for each replicate. Our pipeline does not determine read endedness from a BAM file. You need to explicitly define read endedness. These BAMs should not have unmapped reads or duplicates.
    * Example of 2 singled-ended replicates.
        ```javascript
        {
            "chip.paired_end" : false,
            "chip.nodup_bams" : ["rep1.nodup.bam", "rep2.nodup.bam"]
        }
        ```

4) Starting from TAG-ALIGN BEDs
    * Define a TAG-ALIGN for each replicate. Our pipeline does not determine read endedness from a TAG-ALIGN file. You need to explicitly define read endedness.
    * Example of 4 paired-ended replicates.

        ```javascript
        {
            "chip.paired_end" : true,
            "chip.tas" : ["rep1.tagAlign.gz", "rep2.tagAlign.gz", "rep3.tagAlign.gz", "rep3.tagAlign.gz"]
        }
        ```

You need to define controls for TF ChIP-seq pipeline. Skip this if you want to run histone ChIP-seq pipelines. You can define controls similarly to experiment IP replicates. Just add `ctl_` prefix to parameter names.

1) Control FASTQs
    * Technical replicates for each bio-rep will be **MERGED** in the very early stage of the pipeline. Each read end R1 and R2 have separate arrays `chip.ctl_fastqs_repX_R1` and `chip.ctl_fastqs_repX_R2`. Do not define R2 array for single-ended replicates.
    * Example of 3 paired-ended biological replicates and 2 technical replicates for each bio rep. Two technical replicates `BIOREPX_TECHREP1.R1.fq.gz` and `BIOREPX_TECHREP2.R1.fq.gz` for each bio replicate will be merged.

        ```javascript
        {
            "chip.ctl_paired_end" : true,
            "chip.ctl_fastqs_rep1_R1" : ["BIOREP1_TECHREP1.R1.fq.gz", "BIOREP1_TECHREP2.R1.fq.gz"],
            "chip.ctl_fastqs_rep1_R2" : ["BIOREP1_TECHREP1.R2.fq.gz", "BIOREP1_TECHREP2.R2.fq.gz"],
            "chip.ctl_fastqs_rep2_R1" : ["BIOREP2_TECHREP1.R1.fq.gz", "BIOREP2_TECHREP2.R1.fq.gz"],
            "chip.ctl_fastqs_rep2_R2" : ["BIOREP2_TECHREP1.R2.fq.gz", "BIOREP2_TECHREP2.R2.fq.gz"],
        }
        ```

2) Control BAMs
    * Define a BAM for each replicate. Our pipeline does not determine read endedness from a BAM file. You need to explicitly define read endedness.
    * Example of 3 singled-ended replicates.

        ```javascript
        {
            "chip.ctl_paired_end" : false,
            "chip.ctl_bams" : ["ctl1.bam", "ctl2.bam", "ctl3.bam"]
        }
        ```

3) Control BAMs
    * Define a filtered/deduped BAM for each replicate. Our pipeline does not determine read endedness from a BAM file. You need to explicitly define read endedness. These BAMs should not have unmapped reads or duplicates.
    * Example of 2 singled-ended replicates.
        ```javascript
        {
            "chip.ctl_paired_end" : false,
            "chip.ctl_nodup_bams" : ["ctl1.nodup.bam", "ctl2.nodup.bam"]
        }
        ```

4) Control TAG-ALIGN BEDs
    * Define a TAG-ALIGN for each replicate. Our pipeline does not determine read endedness from a TAG-ALIGN file. You need to explicitly define read endedness.
    * Example of 4 paired-ended replicates.

        ```javascript
        {
            "chip.ctl_paired_end" : true,
            "chip.ctl_tas" : ["ctl1.tagAlign.gz", "ctl2.tagAlign.gz", "ctl3.tagAlign.gz", "ctl4.tagAlign.gz"]
        }
        ```

You can also mix up different data types for individual bio replicate and control. For example, pipeline can start from FASTQs for rep1 (SE) and rep3 (PE), BAMs for rep2 (SE), NODUP_BAMs for rep4 (SE) and TAG-ALIGNs for rep5 (PE). This example has two controls (ctl1: SE BAM, ctl2: PE FASTQs).

```javascript
{
    "chip.paired_ends" : [false, false, true, false, true],
    "chip.fastqs_rep1_R1" : ["rep1.fastq.gz"],
    "chip.fastqs_rep3_R1" : ["rep3.R1.fastq.gz"],
    "chip.fastqs_rep3_R2" : ["rep3.R2.fastq.gz"],
    "chip.bams" : [null, "rep2.bam", null, null, null],
    "chip.nodup_bams" : [null, null, null, "rep4.nodup.bam", null],
    "chip.tas" : [null, null, null, null, "rep5.tagAlign.gz"],

    "chip.ctl_paired_ends" : [false, true],
    "chip.ctl_fastqs_rep2_R1" : ["ctl2.R1.fastq.gz"],    
    "chip.ctl_fastqs_rep2_R2" : ["ctl2.R2.fastq.gz"],    
    "chip.ctl_bams" : ["ctl1.bam", null],
}
```

## Resources

> **WARNING**: It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task. The following parameters are provided for users who want to run our pipeline with Caper's `local` on HPCs and 2).

Resources defined here are **PER BIO REPLICATE**. Therefore, total number of cores will be approximately `chip.align_cpu` x `NUMBER_OF_BIO_REPLICATES` because `align` is a bottlenecking task of the pipeline. This total number of cores will be useful **ONLY** when you use a `local` backend of Caper and manually `qsub` or `sbatch` your job. `disk_factor` is used for GCP/AWS/DNAnexus only.

See [this](input.md#resource-parameters) for details.
