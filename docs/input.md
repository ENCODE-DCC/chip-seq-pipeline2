Input JSON
==========

An input JSON file includes all input parameters and metadata for running pipelines:

1) Reference genome (hg38, mm10, hg19, ...) and genome specific parameters (indices, ...).
2) Input data file paths/URIs (FASTQs, BAMs, TAG-ALIGNs, ...).
3) Pipeline parameters.
4) Resource for instances/jobs.

## Reference genome

We currently support 4 genomes. You can also [build a genome database for your own genome](build_genome_database.md).

|genome|source|built from|
|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|

Choose one TSV file for `"chip.genome_tsv"` in your input JSON. `[GENOME]` should be `hg38`, `mm10`, `hg19` or `mm9`.

|platform|path/URI|
|-|-|
|Google Cloud Platform|`gs://encode-pipeline-genome-data/[GENOME]_google.tsv`|
|DNANexus (CLI)|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:data/pipeline-genome-data/[GENOME]_dx.tsv`|
|DNANexus (CLI, Azure)|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:data/pipeline-genome-data/[GENOME]_dx.tsv`|
|DNANExus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data)|
|DNANExus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data)|
|Stanford Sherlock|`genome/scg/[GENOME]_scg.tsv`|
|Stanford SCG|`genome/sherlock/[GENOME]_sherlock.tsv`|
|Local/SLURM/SGE|You need to [build a genome database](build_genome_database.md). |

## Input data file (IP experiment data)

Choose any data type (FASTQ, BAM, nodup/filtered BAM, TAG-ALIGN and PEAK) you want and DO NOT define arrays for other types. For FASTQs we provide two ways to define them since DNANexus web UI supports up to an 1-dim array. Choose between 3-dim `fastqs` or 1-dim `fastqs_rep[REP_ID]_R[READ_END_ID]` according to your preference. The pipeline supports up to 6 replicates.

* `"chip.fastqs"` : 3-dimensional array with FASTQ file path/URI.
    - 1st dimension: replicate ID
    - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
    - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
* `"chip.fastqs_rep1_R1"` : Array of FASTQ file to be merged for rep1-R1.
* `"chip.fastqs_rep1_R2"` : Array of FASTQ file to be merged for rep1-R2. Do not define if your FASTQ is single ended.
* `"chip.fastqs_rep2_R1"` : Array of FASTQ file to be merged for rep2-R1. Do not define if you don't have replicate 2.
* `"chip.fastqs_rep2_R2"` : Array of FASTQ file to be merged for rep2-R2. Do not define if you don't have replicate 2.
* `"chip.fastqs_rep3_R1"` : Array of FASTQ file to be merged for rep3-R1. Do not define if you don't have replicate 3.
* `"chip.fastqs_rep3_R2"` : Array of FASTQ file to be merged for rep3-R2. Do not define if you don't have replicate 3.
* `"chip.fastqs_rep4_R1"` : Array of FASTQ file to be merged for rep4-R1. Do not define if you don't have replicate 4.
* `"chip.fastqs_rep4_R2"` : Array of FASTQ file to be merged for rep4-R2. Do not define if you don't have replicate 4.
* `"chip.bams"` : Array of raw (unfiltered) BAM file path/URI.
    - 1st dimension: replicate ID
* `"chip.nodup_bams"` : Array of filtered (deduped) BAM file path/URI.
    - 1st dimension: replicate ID
* `"chip.tas"` : Array of TAG-ALIGN file path/URI.
    - 1st dimension: replicate ID
* `"chip.peaks"` : Array of NARROWPEAK file path/URI.
    - 1st dimension: replicate ID
* `"chip.peaks_pr1"` : Array of NARROWPEAK file path/URI for 1st self pseudo replicate of replicate ID.
    - 1st dimension: replicate ID
* `"chip.peaks_pr2"` : Array of NARROWPEAK file path/URI for 2nd self pseudo replicate of replicate ID.
    - 1st dimension: replicate ID
* `"chip.peak_ppr1"` : NARROWPEAK file path/URI for pooled 1st pseudo replicates.
* `"chip.peak_ppr2"` : NARROWPEAK file path/URI for pooled 2nd pseudo replicates.
* `"chip.peak_pooled"` : NARROWPEAK file path/URI for pooled replicate.

If starting from peaks then always define `"chip.peaks"`. Define `"chip.peaks_pr1"`, `"chip.peaks_pr2"`, `"chip.peak_pooled"`, `"chip.peak_ppr1"` and `"chip.peak_ppr2"` according to the following rules:

```
if num_rep>1:
    if true_rep_only: peak_pooled, 
    else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
else:
    if true_rep_only: "not the case!"
    else: peaks_pr1[], peaks_pr2[]
```

## Input data file (Control data)

* `"chip.ctl_fastqs"` : 3-dimensional array with FASTQ file path/URI.
    - 1st dimension: replicate ID
    - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
    - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
* `"chip.ctl_fastqs_rep1_R1"` : Array of FASTQ file to be merged for rep1-R1.
* `"chip.ctl_fastqs_rep1_R2"` : Array of FASTQ file to be merged for rep1-R2. Do not define if your FASTQ is single ended.
* `"chip.ctl_fastqs_rep2_R1"` : Array of FASTQ file to be merged for rep2-R1. Do not define if you don't have replicate 2.
* `"chip.ctl_fastqs_rep2_R2"` : Array of FASTQ file to be merged for rep2-R2. Do not define if you don't have replicate 2.
* `"chip.ctl_fastqs_rep3_R1"` : Array of FASTQ file to be merged for rep3-R1. Do not define if you don't have replicate 3.
* `"chip.ctl_fastqs_rep3_R2"` : Array of FASTQ file to be merged for rep3-R2. Do not define if you don't have replicate 3.
* `"chip.ctl_fastqs_rep4_R1"` : Array of FASTQ file to be merged for rep4-R1. Do not define if you don't have replicate 4.
* `"chip.ctl_fastqs_rep4_R2"` : Array of FASTQ file to be merged for rep4-R2. Do not define if you don't have replicate 4.
* `"chip.ctl_bams"` : Array of raw (unfiltered) BAM file path/URI.
    - 1st dimension: replicate ID
* `"chip.ctl_nodup_bams"` : Array of filtered (deduped) BAM file path/URI.
    - 1st dimension: replicate ID
* `"chip.ctl_tas"` : Array of TAG-ALIGN file path/URI.
    - 1st dimension: replicate ID

## Pipeline parameters

1. General

    Choose pipeline type: TF (`tf`) or Histone (`histone`) ChIP-Seq. Default peak caller for TF and Histone ChIP-Seq pipelines are `spp` and `macs2`, respectively. However you can also manually specify a peak caller for these pipeline types. MACS2 can work without controls but SPP cannot. Therefore, if a peak caller is chosen as `spp` then make sure to define control data.

    * `"chip.pipeline_type` : `tf` for TF ChIP-Seq. `histone` for Histone ChIP-Seq.
    * `"chip.peak_caller` (optional) : Choose between `macs2` or `spp` if you don't want to use a default peak caller for `pipeline_type` chosen.

    Input data endedness.

    * `"chip.paired_end"` : Set it as `true` if input data are paired end, otherwise `false`.

    Other optional settings.

    * `"chip.align_only"` : (optional) Disable all downstream analysis (peak calling, ...) after mapping.
    * `"chip.true_rep_only"` : (optional) Set it as `true` to disable all analyses (including IDR, naive-overlap and reproducibility QC) related to pseudo replicates. This flag suppresses `"chip.enable_idr"`.
    * `"chip.disable_xcor` : (optional) Disable cross-correlation analysis.

    * `"chip.title"` : (optional) Name of sample.
    * `"chip.description"` : (optional) Description for sample.

2. Trim FASTQ settings (for paired end dataset only).

    * `"chip.trim_bp"` : For paired end dataset only. Number of basepairs after trimming FASTQ. It’s 50 by default. Trimmed FASTQS is only used for cross-correlation analysis. FASTQ mapping is not affected by this parameter.

3. Filter/dedup (post-alignment) settings (remove a prefix `chip.` for DNANexus CLI).

    * `"chip.dup_marker"` : (optional) Dup marker. Choose between `picard` (default) and `sambamba`.
    * `"chip.mapq_thresh"` : (optional) Threshold for low MAPQ reads removal (default: 30).
    * `"chip.no_dup_removal"` : (optional) No dup reads removal when filtering BAM.

4. BAM-2-TAGALIGN settings (remove a prefix `chip.` for DNANexus CLI).

    Pipeline filters out chrM reads by default.

    * `"chip.regex_filter_reads"` : (optional) Perl-style regular expression pattern to remove matching reads from TAGALIGN (default: `chrM`).
    * `"chip.subsample_reads"` : (optional) Number of reads to subsample TAGALIGN. Subsampled TAGALIGN will be used for all downstream analysis (MACS2, IDR, naive-overlap).

5. Choose control settings.

    * `"chip.ctl_depth_ratio"` : (optional) if ratio between controls is higher than this then always use pooled control for all exp rep (default: 1.2).
    * `"chip.always_use_pooled_ctl"` : (optional) Always use pooled control for all exp replicates (ignoring ctl_depth_ratio).

6. Cross correlation analysis settings (remove a prefix `chip.` for DNANexus CLI).

    For paired end FASTQ data set, only one read end (R1) will be trimmed and used for this analysis.

    * `"chip.xcor_subsample_reads"` : (optional) Number of reads to subsample TAGALIGN. This will not affect downstream analysis.

7. MACS2 settings

    **DO NOT DEFINE MACS2 PARAMETERS IN `"chip.macs2"` SCOPE**. All MACS2 parameters must be defined in `"chip"` scope.

    * `"chip.macs2_cap_num_peak"` : (optional) Cap number of raw peaks called from MACS2 (default: 500000).
    * `"chip.pval_thresh"` : (optional) P-value threshold (default: 0.01).

8. SPP settings

    **DO NOT DEFINE SPP PARAMETERS IN `"chip.spp"` SCOPE**. All SPP parameters must be defined in `"chip"` scope.

    * `"chip.spp_cap_num_peak"` : (optional) Cap number of raw peaks called from SPP (default: 300000).

9. IDR settings

    **DO NOT DEFINE IDR PARAMETERS IN `"chip.idr"` SCOPE**. All IDR parameters must be defined in `"chip"` scope.

    * `"chip.enable_idr"` : (optional) Set it as `true` to enable IDR on raw peaks.
    * `"chip.idr_thresh"` : (optional) IDR threshold (default: 0.05).

## Resource

**RESOURCES DEFINED IN AN INPUT JSON ARE PER TASK**. For example, if you have FASTQs for 2 replicates (2 tasks) and set `cpu` for `bwa` task as 4 then total number of cpu cores to map FASTQs is 2 x 4 = 8.

CPU (`cpu`), memory (`mem_mb`) settings are used for submitting jobs to cluster engines (SGE and SLURM) and Cloud platforms (Google Cloud Platform, AWS, ...). VM instance type on cloud platforms will be automatically chosen according to each task's `cpu` and `mem_mb`. Number of cores for tasks without `cpu` parameter is fixed at 1.

* `"chip.bwa_cpu"` : (optional) Number of cores for `bwa` (default: 4).
* `"chip.filter_cpu"` : (optional) Number of cores for `filter` (default: 2).
* `"chip.bam2ta_cpu"` : (optional) Number of cores for `bam2ta` (default: 2).
* `"chip.xcor_cpu"` : (optional) Number of cores for `xcor` (default: 2).
* `"chip.trim_adapter_mem_mb"` : (optional) Max. memory limit in MB for `trim_adapter` (default: 10000).
* `"chip.bwa_mem_mb"` : (optional) Max. memory limit in MB for `bwa` (default: 20000).
* `"chip.filter_mem_mb"` : (optional) Max. memory limit in MB for `filter` (default: 20000).
* `"chip.bam2ta_mem_mb"` : (optional) Max. memory limit in MB for `bam2ta` (default: 10000).
* `"chip.spr_mem_mb"` : (optional) Max. memory limit in MB for `spr` (default: 12000).
* `"chip.xcor_mem_mb"` : (optional) Max. memory limit in MB for `xcor` (default: 10000).
* `"chip.macs2_mem_mb"` : (optional) Max. memory limit in MB for `macs2` (default: 16000).

Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

* `"chip.bwa_disks"` : (optional) Disks for `bwa` (default: "local-disk 100 HDD").
* `"chip.filter_disks"` : (optional) Disks for `filter` (default: "local-disk 100 HDD").
* `"chip.bam2ta_disks"` : (optional) Disks for `bam2ta` (default: "local-disk 100 HDD").
* `"chip.xcor_disks"` : (optional) Disks for `xcor` (default: "local-disk 100 HDD").
* `"chip.macs2_disks"` : (optional) Disks for `macs2` (default: "local-disk 100 HDD").

Walltime (`time`) settings (for SGE and SLURM only).

* `"chip.bwa_time_hr"` : (optional) Walltime for `bwa` (default: 48).
* `"chip.filter_time_hr"` : (optional) Walltime for `filter` (default: 24).
* `"chip.bam2ta_time_hr"` : (optional) Walltime for `bam2ta` (default: 6).
* `"chip.xcor_time_hr"` : (optional) Walltime for `xcor` (default: 24).
* `"chip.macs2_time_hr"` : (optional) Walltime for `macs2` (default: 24).

