# Input JSON

An input JSON file is a file which must include all the information needed to run this pipeline. Hence, it must include the absolute paths to all the control and experimental fastq files; paths to all the genomic data files needed for this pipeline, and it must also specify the parameters and the metadata needed for running this pipeline. If the parameters are not specified in an input JSON file, default values will be used. We provide a set of template JSON files: [minimum](../example_input_json/template.json) and [full](../example_input_json/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

>**IMPORTANT**: ALWAYS USE ABSOLUTE PATHS.

# Checklist

Mandatory parameters.

1) Pipeline type
    * `chip.pipeline_type`: `tf` for TF ChIP-seq or `histone` for histone ChIP-seq. One major difference between two types is that `tf` uses `spp` peak caller with controls but `histone` uses `macs2` peak caller without controls. 

2) Experiment title/description
    * `chip.title`: experiment title for a final HTML report.
    * `chip.description`: experiment description for a final HTML report.

3) Read endedness
    * `chip.paired_end`: `true` if ALL replicates are paired-ended.
    * (Optional) `chip.paired_ends`: For samples with mixed read ends, you can define read endedness for each biological replicate (e.g. `[true, false]` means paired-ended biorep-1 and single-ended biorep-2).
    * `chip.ctl_paired_end`: `true` if ALL controls are paired-ended. If not defined then `chip.paired_end` will be used.
    * (Optional) `chip.ctl_paired_ends`: For controls with mixed read ends, you can define read endedness for each biological replicate (e.g. `[true, false]` means paired-ended biorep-1 and single-ended biorep-2). If not defined then `chip.paired_ends` will be used.

4) Reference genome
    * `chip.genome_tsv`: Use `https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/[GENOME]_caper.tsv`.
    * Supported `GENOME`s: are hg38, mm10, hg19 and mm9.
    * We provide a genome TSV file that defines all genome-specific parameters and reference data files. Caper will automatically download big reference data files from our ENCODE repository.
    * However, we also have reference data mirrors for [some platforms](input_details.md/#reference-genome) (GCP, AWS, Sherlock, SCG, ...). On these platforms, you can use a different TSV file to prevent downloading such big reference data.
    * To build a new TSV file from use your own FASTA (`.fa` and `.2bit`) see [this](build_genome_database.md).

5) [Input files](#input-files) and [adapters](#adapters)
    * See [this](#input-files) for how to define FASTQ/BAM/TAG-ALIGNs for your sample.
    * See [this](#adapters) for how to define adapters to be trimmed.

6) Important parameters
    * `chip.crop_length`: Crop FASTQs with Trimmomatic. **WARNING**: Check your FASTQs' read length first. Reads SMALLER than this will be excluded from mapping, hence not included in output BAM files and all downstream analyses. It is 0 (disabled) by default.
    * `chip.always_use_pooled_ctl`: (For TF ChIP-seq only) Always use a pooled control to compare with each replicate. If a single control is given then use it. It is disabled by default.
    * `chip.ctl_depth_ratio`: (For TF ChIP-seq only) If ratio of depth between controls is higher than this. then always use a pooled control for all replicates. It's 1.2 by default.

7) [Resources](#resources)
    * If your FASTQs/BAMs are big (>10GB) then try with higher resource settings, especially for memory (`chip.[TASK_NAME]_mem_mb`).

Optional parameters.

8) Useful parameters
    * `chip.subsample_reads`: Subsample experimet reads. This will affect all downsteam analyses including peak-calling. It's 0 by default, which means no subsampling.
    * `chip.ctl_subsample_reads`: Subsample control reads. This will affect all downsteam analyses including peak-calling. It's 0 by default, which means no subsampling.
    * `chip.fraglen`: Array of Integers. Fragment length for each bio replicate. If you start from FASTQs then our pipeline automatically estimate it from cross-correlation analysis (task `xcor`) result since such analysis requires a special treamtment for FASTQs. It is possible that fragment length is not estimated correctly (or pipeline can fail due to negative fraglen) if you start from different types (BAM/TAG-ALIGN). For such case, you can manually define fragment length for each bio rep. (e.g. `[200, 150]` means 200 for rep1 and 150 for rep2).

9) Flags
    * `chip.align_only`: Peak calling and its downstream analyses will be disabled. Useful if you just want to align your FASTQs into filtered BAMs/TAG-ALIGNs and don't want to call peaks on them.
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

Resources defined here are **PER BIO REPLICATE**. Therefore, total number of cores will be approximately `chip.align_cpu` x `NUMBER_OF_BIO_REPLICATES` because `align` is a bottlenecking task of the pipeline. This total number of cores will be useful **ONLY** when you use a `local` backend of Caper and manually `qsub` or `sbatch` your job. `disks` is used for Google Cloud and DNAnexus only.


Parameter|Default
---------|-------
`chip.align_cpu` | 4
`chip.align_mem_mb` | 20000
`chip.align_time_hr` | 48
`chip.align_disks` | `local-disk 400 HDD`

Parameter|Default
---------|-------
`chip.filter_cpu` | 2
`chip.filter_mem_mb` | 20000
`chip.filter_time_hr` | 24
`chip.filter_disks` | `local-disk 400 HDD`

Parameter|Default
---------|-------
`chip.bam2ta_cpu` | 2
`chip.bam2ta_mem_mb` | 10000
`chip.bam2ta_time_hr` | 6
`chip.bam2ta_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`chip.spr_mem_mb` | 16000

Parameter|Default
---------|-------
`chip.jsd_cpu` | 2
`chip.jsd_mem_mb` | 12000
`chip.jsd_time_hr` | 6
`chip.jsd_disks` | `local-disk 200 HDD`

Parameter|Default
---------|-------
`chip.xcor_cpu` | 2
`chip.xcor_mem_mb` | 16000
`chip.xcor_time_hr` | 24
`chip.xcor_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`chip.call_peak_cpu` | 2
`chip.call_peak_mem_mb` | 16000
`chip.call_peak_time_hr` | 24
`chip.call_peak_disks` | `local-disk 200 HDD`

Parameter|Default
---------|-------
`chip.macs2_signal_track_mem_mb` | 16000
`chip.macs2_signal_track_time_hr` | 24
`chip.macs2_signal_track_disks` | `local-disk 400 HDD`

> **IMPORTANT**: If you see Java memory errors, check the following resource parameters.

There are special parameters to control maximum Java heap memory (e.g. `java -Xmx4G`) for Picard tools. They are strings including size units. Such string will be directly appended to Java's parameter `-Xmx`.

Parameter|Default
---------|-------
`chip.filter_picard_java_heap` | = `chip.filter_mem_mb`
`chip.align_trimmomatic_java_heap` | = `chip.align_mem_mb`
`chip.gc_bias_picard_java_heap` | `10G`
