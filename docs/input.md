# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file. We provide a set of template JSON files: [minimum](../example_input_json/template.json) and [full](../example_input_json/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

Please read through the following step-by-step instruction to compose a input JSON file.

>**IMPORTANT**: ALWAYS USE ABSOLUTE PATHS.

## Pipeline metadata

Parameter|Description
---------|-----------
`chip.title`| Title for experiment which will be shown in a final HTML report
`chip.description`| Description for experiment which will be shown in a final HTML report

## Pipeline parameters

Parameter|Default|Description
---------|-------|-----------
`chip.pipeline_type`| `tf` | `tf` for TF ChIP-seq, `histone` for Histone ChIP-seq and `control` for mapping control FASTQs only. If it is `control` then `chip.align_only` is automatically turned on and pipeline will align FASTQs defined in `chip.fastqs_repX_RY` (where `X` means control replicate ID and `Y` is read-end) as controls. For control mode, do not define  `chip.ctl_fastqs_repX_RY`. Pipeline will skip cross-correlation/JSD/GC-bias analysis for control mode.
`chip.align_only`| false | Peak calling and its downstream analyses will be disabled. Useful if you just want to map your FASTQs into filtered BAMs/TAG-ALIGNs and don't want to call peaks on them. Even though `chip.pipeline_type` does not matter for align only mode, you still need to define it since it is a required parameter in WDL. Define it as `tf` for such cases.
`chip.true_rep_only` | false | Disable pseudo replicate generation and all related analyses

## Reference genome

All reference genome specific reference files/parameters can be defined in a single TSV file `chip.genome_tsv`. However, you can also individally define each file/parameter instead of a TSV file. If both a TSV file and individual parameters are defined, then individual parameters will override those defined in a TSV file. For example, if you define both `chip.genome_tsv` and `chip.blacklist`, then `chip.blacklist` will override that is defined in `chip.genome_tsv`. This is useful when you want to use your own for a specific parameter while keeping all the other parameters same as original.

Parameter|Type|Description
---------|-------|-----------
`chip.genome_tsv`| File | Choose one of the TSV files listed below or build your own
`chip.genome_name`| String | Name of genome (e.g. hg38, hg19, ...)
`chip.ref_fa`| File | Reference FASTA file
`chip.ref_mito_fa`| File | Mito-only reference FASTA file
`chip.bowtie2_idx_tar`| File | Bowtie2 index TAR file (uncompressed) built from FASTA file
`chip.bowtie2_mito_idx_tar`| File | Mito-only Bowtie2 index TAR file (uncompressed) built from FASTA file
`chip.bwa_idx_tar`| File | BWA index TAR file (uncompressed) built from FASTA file with `bwa index`
`chip.bwa_mito_idx_tar`| File | Mito-only BWA index TAR file (uncompressed) built from FASTA file with `bwa index`
`chip.chrsz`| File | 2-col chromosome sizes file built from FASTA file with `faidx`
`chip.blacklist`| File | 3-col BED file. Peaks overlapping these regions will be filtered out
`chip.blacklist2`| File | Second blacklist. Two blacklist files (`atac.blacklist` and `atac.blacklist2`) will be merged.
`chip.gensz`| String | MACS2's genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)
`chip.mito_chr_name`| String | Name of mitochondrial chromosome (e.g. chrM)
`chip.regex_bfilt_peak_chr_name`| String | Perl style reg-ex to keep peaks on selected chromosomes only matching with this pattern (default: `chr[\dXY]+`. This will keep chr1, chr2, ... chrX and chrY in `.bfilt.` peaks file. chrM is not included here)

We currently provide TSV files for 4 genomes as shown in the below table. `GENOME` should be `hg38`, `mm10`, `hg19` or `mm9`. You can [download/build](build_genome_database.md) it on your local computer. You can also [build a genome database for your own genome](build_genome_database.md).

Platform|Path/URI
-|-
Google Cloud Platform|`gs://encode-pipeline-genome-data/genome_tsv/v1/[GENOME]_gcp.tsv`
Stanford Sherlock|`/home/groups/cherry/encode/pipeline_genome_data/genome_tsv/v1/[GENOME]_sherlock.tsv`
Stanford SCG|`/reference/ENCODE/pipeline_genome_data/genome_tsv/v1/[GENOME]_scg.tsv`
Local/SLURM/SGE|You need to [build](build_genome_database.md) or [download]() a genome database]. 
DNAnexus (CLI)|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v1/[GENOME]_dx.tsv`
DNAnexus (CLI, Azure)|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v1/[GENOME]_dx_azure.tsv`
DNAnexus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v1)
DNAnexus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v1)

Additional information about each genome:

|Genome|Source|built from|
|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|

## How to download genome database

1. Choose `GENOME` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory.
    ```bash
    $ bash genome/download_genome_data.sh [GENOME] [DESTINATION_DIR]
    ```
2. Find a TSV file on the destination directory and use it for `"chip.genome_tsv"` in your input JSON.

## Input genomic data

Choose endedness of your dataset first.

Parameter|Description
---------|-----------
`chip.paired_end`| Boolean to define endedness for ALL IP replicates. This will override per-replicate definition in `chip.paired_ends`
`chip.ctl_paired_end`| Boolean to define endedness for ALL control replicates. This will override per-replicate definition in `chip.ctl_paired_ends`
`chip.paired_ends`| Array of Boolean to define endedness for each replicate
`chip.ctl_paired_ends`| Array of Boolean to define endedness for each control replicate

Define `chip.paired_end` and `chip.ctl_paired_end` if all replicates (or control replicates) in your dataset has the same endedness. You can also individually define endedness for each replicate and control replicate. For example, rep1, rep2 are PE and rep3 is SE. control rep1 is SE and control rep2 is PE.

```javascript
{
    "chip.paired_ends" : [true, true, false],
    "chip.ctl_paired_ends" : [false, true]
}
```

Pipeline can start from any of the following data type (FASTQ, BAM, NODUP_BAM and TAG-ALIGN).

Parameter|Description
---------|-----------
`chip.fastqs_repX_R1`| Array of R1 FASTQ files for replicate X. These files will be merged into one FASTQ file for rep X.
`chip.fastqs_repX_R2`| Array of R2 FASTQ files for replicate X. These files will be merged into one FASTQ file for rep X. Do not define for single ended dataset. 
`chip.bams`| Array of BAM file for each replicate. (e.g. `["rep1.bam", "rep2.bam", ...]`)
`chip.nodup_bams`| Array of filtered/deduped BAM file for each replicate.
`chip.tas`| Array of TAG-ALIGN file for each replicate.

For controls:

Parameter|Description
---------|-----------
`chip.ctl_fastqs_repX_R1`| Array of R1 FASTQ files for control replicate X. These files will be merged into one FASTQ file for rep X.
`chip.ctl_fastqs_repX_R2`| Array of R2 FASTQ files for control replicate X. These files will be merged into one FASTQ file for rep X. Do not define for single ended dataset. 
`chip.ctl_bams`| Array of BAM file for each control replicate. (e.g. `["ctl_rep1.bam", "ctl_rep2.bam", ...]`)
`chip.ctl_nodup_bams`| Array of filtered/deduped BAM file for each control replicate.
`chip.ctl_tas`| Array of TAG-ALIGN file for each control replicate.

You can mix up different data types for individual replicate/control replicate. For example, pipeline can start from FASTQs for rep1 and rep3, BAMs for rep2, NODUP_BAMs for rep4 and TAG-ALIGNs for rep5. You can define similarly for control replicates.

```javascript
{
    "chip.fastqs_rep1_R1" : ["rep1.fastq.gz"],
    "chip.fastqs_rep3_R1" : ["rep3.fastq.gz"],
    "chip.bams" : [null, "rep2.bam", null, null, null],
    "chip.nodup_bams" : [null, null, null, "rep4.nodup.bam", null],
    "chip.tas" : [null, null, null, null, "rep5.tagAlign.gz"]
}
```

## Optional mapping parameters

Parameter|Type|Default|Description
---------|---|----|-----------
`chip.aligner` | String | bowtie2 | Currently supported aligners: bwa and bowtie2.
`chip.crop_length` | Int | 0 | Crop FASTQs with Trimmomatic (using parameters `CROP`). 0: cropping disabled.
`chip.crop_length_tol` | Int | 2 | Trimmomatic's `MINLEN` will be set as `chip.crop_length` - `abs(chip.crop_length_tol)` where reads shorter than `MINLEN` will be removed, hence not included in output BAM files and all downstream analyses.
`chip.use_bwa_mem_for_pe` | Boolean | false | For PE dataset, uise bwa mem instead of bwa aln.


## Optional filtering parameters

Parameter|Default|Description
---------|-------|-----------
`chip.mapq_thresh` | 30 for bwa, 255 for bowtie2 | Threshold for mapped reads quality (samtools view -q). If not defined, automatically determined according to aligner.
`chip.dup_marker` | `picard` | Choose a dup marker between `picard` and `sambamba`. `picard` is recommended, use `sambamba` only when picard fails.
`chip.no_dup_removal` | false | Skip dup removal in a BAM filtering stage.

## Optional subsampling parameters

Parameter|Default|Description
---------|-------|-----------
`chip.subsample_reads` | 0 | Subsample reads (0: no subsampling). For PE dataset, this is not a number of read pairs but number of reads. Subsampled reads will be used for all downsteam analyses including peak-calling. Subsampling is done during BAM to TAG-ALIGN conversion. Output TAG-ALIGN and all donwsteam analyses like peak-calling will be affected.
`chip.ctl_subsample_reads` | 0 | Subsample control reads (not number of read pairs). For PE dataset, this is not a number of read pairs but number of reads. Subsampling is done during BAM to TAG-ALIGN conversion. Output control TAG-ALIGN and all donwsteam analyses like peak-calling will be affected.
`chip.xcor_subsample_reads` | 15000000 | Subsample reads for cross-corr. analysis only (0: no subsampling). Subsampled reads will be used for cross-corr. analysis only
`chip.ctl_depth_limit` | 200000000 | Hard limit for automatic subsampling control. For PE dataset, this is not a number of read pairs but number of reads. This is different from manual control subsampling controlled by `chip.xcor_subsample_reads`. This is an additional subsampling for controls in a peak-calling task.
`chip.exp_ctl_depth_ratio_limit` | 5.0 | For each experiment replicate, corresponding chosen control will be subsampled if control's read is deeper than replicate's multiplied by this number.


## Optional cross-correlation analysis parameters

Parameter|Default|Description
---------|-------|-----------
`chip.xcor_trim_bp` | 50 | Trim R1 fastq for cross-correlation analysis only. Trimmed fastqs will not be used for any other analyses
`chip.use_filt_pe_ta_for_xcor` | false | Use filtered PE BAM/TAG-ALIGN for cross-correlation analysis ignoring the above trimmed R1 fastq
`chip.xcor_exclusion_range_min` | -500 | Exclusion minimum for cross-corr. analysis. See [description for `-x=<min>:<max>`](https://github.com/kundajelab/phantompeakqualtools) for details. Make sure that it's consistent with default strand shift `-s=-500:5:1500` in `phantompeakqualtools`.
`chip.xcor_exclusion_range_max` |  | Exclusion minimum for cross-corr. analysis. If not defined default value of `max(read length + 10, 50)` for TF and `max(read_len + 10, 100)` for histone are used.

## Optional control parameters

Parameter|Default|Description
---------|-------|-----------
`chip.always_use_pooled_ctl` | true | Choosing an appropriate control for each replicate. Always use a pooled control to compare with each replicate. If a single control is given then use it.
`chip.ctl_depth_ratio` | 1.2 | If ratio of depth between controls is higher than this. then always use a pooled control for all replicates.

## Optional peak-calling parameters

Parameter|Default|Description
---------|-------|-----------
`chip.peak_caller`| `spp` for `tf` type<br>`macs2` for `histone` type| `spp` or `macs2`. `spp` requires control<br>`macs2` can work without controls
`chip.cap_num_peak` | 500000 (`macs2`), 300000 (`spp`) | Cap number of peaks called.
`chip.pval_thresh` | 0.01 | P-value threshold for peak-caller MACS2 (macs2 callpeak -p).
`chip.idr_thresh` | 0.05 | IDR (irreproducible discovery rate) threshold.
`chip.fdr_thresh` | 0.01 | FDR threshold for peak-caller SPP (run_spp.R -fdr=).

## Optional pipeline flags

Parameter|Default|Description
---------|-------|-----------
`chip.enable_jsd` | true | Enable deeptools fingerprint (JS distance)
`chip.enable_gc_bias` | true | Enable GC bias calculation
`chip.enable_count_signal_track` | false | Enable count signal track generation

## Optional parameter for fragment length

Our pipeline automatically estimate fragment lengths (required for TF ChIP-Seq) from cross-correlation (task `xcor`) anaylses, but `chip.fraglen` will override those estimated ones. Use this if your pipeline fails due to invalid (negative) fragment length estimated from the cross-correlation analysis.

Parameter|Type | Description
---------|-----|-----------
`chip.fraglen` | `Array[Int]` | Fragment length for each replicate.

## Other optional parameters

Parameter|Default|Description
---------|-------|-----------
`chip.filter_chrs` | `[]` (empty array of string) | Array of chromosome names to be filtered out from a final (filtered/nodup) BAM. No chromosomes are filtered out by default.

## Resource parameters

> **WARNING**: It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task. The following parameters are provided for users who want to run our pipeline with Caper's `local` on HPCs and 2).

Resources defined here are PER REPLICATE. Therefore, total number of cores will be MAX(`chip.align_cpu` x `NUMBER_OF_REPLICATES`, `chip.call_peak_cpu` x 2 x `NUMBER_OF_REPLICATES`) because `align` and `call_peak` (especially for `spp`) are bottlenecking tasks of the pipeline. Use this total number of cores if you manually `qsub` or `sbatch` your job (using local mode of Caper). `disks` is used for Google Cloud and DNAnexus only.

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

There are special parameters to control maximum Java heap memory (e.g. `java -Xmx4G`) for Picard tools. They are strings including size units. Such string will be directly appended to Java's parameter `-Xmx`. If these parameters are not defined then pipeline uses 90% of each task's memory (e.g. `chip.filter_mem_mb`).

Parameter|Default
---------|-------
`chip.filter_picard_java_heap` | = `chip.filter_mem_mb`
`chip.align_trimmomatic_java_heap` | = `chip.align_mem_mb`
`chip.gc_bias_picard_java_heap` | `10G`
