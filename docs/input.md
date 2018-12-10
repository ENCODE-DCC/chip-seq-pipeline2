# Input JSON

An input JSON file includes all input parameters and metadata for running pipelines. Items 1) and 2) are mandatory. Items 3) and 4) are optional so that our pipeline will use default values if they are not defined. However, 

* Mandatory

1. Reference genome.
2. Input data file paths/URIs.

* Optional

3. Pipeline parameters.
4. Resource settings for jobs.

## Templates

We provide two template JSON files for both single ended and paired-end samples. We recommend to use one of these input JSON files instead of that used in the tutorial section. These template JSON files include all parameters of the pipeline with default values defined.

* [template](../examples/template_se.json) for single ended sample
* [template](../examples/template_pe.json) for paired-end sample

Let us take a close look at the following template JSON. Comments are not allowed in a JSON file but we added some comments to help you understand each parameter.
```javascript
{
    ////////// 1) Reference genome //////////
    // Stanford servers: [GENOME]=hg38,hg19,mm10,mm9
    //   Sherlock: /home/groups/cherry/encode/pipeline_genome_data/[GENOME]_sherlock.tsv
    //   SCG4: /reference/ENCODE/pipeline_genome_data/[GENOME]_scg.tsv

    // Cloud platforms (Google Cloud, DNANexus): [GENOME]=hg38,hg19,mm10,mm9
    //   Google Cloud: gs://encode-pipeline-genome-data/[GENOME]_google.tsv
    //   DNANexus: dx://project-BKpvFg00VBPV975PgJ6Q03v6:data/pipeline-genome-data/[GENOME]_dx.tsv
    //   DNANexus(Azure): dx://project-F6K911Q9xyfgJ36JFzv03Z5J:data/pipeline-genome-data/[GENOME]_dx_azure.tsv

    // On other computers download or build reference genome database and pick a TSV from [DEST_DIR].
    //   Downloader: ./genome/download_genome_data.sh [GENOME] [DEST_DIR]
    //   Builder (Conda required): ./conda/build_genome_data.sh [GENOME] [DEST_DIR]

    "chip.genome_tsv" : "/path_to_genome_data/hg38/hg38.tsv",

    ////////// 2) Input data files paths/URIs //////////

    // Read endedness
    "chip.paired_end" : true,

    // If you start from FASTQs then define these, otherwise remove from this file.
    // You can define up to 6 replicates.
    // FASTQs in an array will be merged after trimming adapters.
    // For example, 
    // "rep1_R1_L1.fastq.gz", "rep1_R1_L2.fastq.gz" and "rep1_R1_L3.fastq.gz" will be merged together.
    "chip.fastqs_rep1_R1" : [ "rep1_R1_L1.fastq.gz", "rep1_R1_L2.fastq.gz", "rep1_R1_L3.fastq.gz" ],
    "chip.fastqs_rep1_R2" : [ "rep1_R2_L1.fastq.gz", "rep1_R2_L2.fastq.gz", "rep1_R2_L3.fastq.gz" ],
    "chip.fastqs_rep2_R1" : [ "rep2_R1_L1.fastq.gz", "rep2_R1_L2.fastq.gz" ],
    "chip.fastqs_rep2_R2" : [ "rep2_R2_L1.fastq.gz", "rep2_R2_L2.fastq.gz" ],

    // Define if you have control FASTQs otherwise remove from this file.
    "chip.ctl_fastqs_rep1_R1" : [ "ctl1_R1.fastq.gz" ],
    "chip.ctl_fastqs_rep1_R2" : [ "ctl1_R2.fastq.gz" ],
    "chip.ctl_fastqs_rep2_R1" : [ "ctl2_R1.fastq.gz" ],
    "chip.ctl_fastqs_rep2_R2" : [ "ctl2_R2.fastq.gz" ],

    // If you start from BAMs then define these, otherwise remove from this file.
    // You can define up to 6 replicates. The following example array has two replicates.
    "chip.bams" : [
        "raw_rep1.bam",
        "raw_rep2.bam"
    ],
    // Define if you have control BAMs otherwise remove from this file.
    "chip.ctl_bams" : [
        "raw_ctl1.bam",
        "raw_ctl2.bam"
    ],

    // If you start from filtered/deduped BAMs then define these, otherwise remove from this file.
    // You can define up to 6 replicates. The following example array has two replicates.
    "chip.nodup_bams" : [
        "nodup_rep1.bam",
        "nodup_rep2.bam"
    ],
    // Define if you have control filtered/deduped BAMs otherwise remove from this file.
    "chip.ctl_nodup_bams" : [
        "nodup_ctl1.bam",
        "nodup_ctl2.bam"
    ],

    // If you start from TAG-ALIGNs then define these, otherwise remove from this file.
    // You can define up to 6 replicates. The following example array has two replicates.
    "chip.tas" : [
        "rep1.tagAlign.gz",
        "rep2.tagAlign.gz"
    ],
    // Define if you have control TAG-ALIGNs otherwise remove from this file.
    "chip.ctl_tas" : [
        "ctl1.tagAlign.gz",
        "ctl2.tagAlign.gz"
    ],

    ////////// 3) Pipeline parameters //////////

    // Pipeline title and description
    "chip.title" : "Example (single-ended)",
    "chip.description" : "This is an template input JSON for single-ended sample.",

    // Pipeline type (tf or histone).
    // default peak_caller: spp for tf, macs2 for histone
    "chip.pipeline_type" : "tf",
    // You can also manually specify a peak_caller
    "chip.peak_caller" : "spp",

    // Pipeline will not proceed to post alignment steps (peak-calling, ...).
    // You will get QC report for alignment only.
    "chip.align_only" : false,
    "chip.true_rep_only" : false,

    // Disable deeptools fingerprint (JS distance)
    "chip.disable_fingerprint" : false,

    // Trim paired ended fastqs for cross-correlation analysis only
    // Trimmed fastqs will not be used for any other analyses
    "chip.trim_bp" : 50,

    // Choose a dup marker between picard and sambamba
    // picard is recommended, use sambamba only when picard fails.
    "chip.dup_marker" : "picard",

    // Threshold for mapped reads quality (samtools view -q)
    "chip.mapq_thresh" : 30,

    // Skip dup removal in a BAM filtering stage.
    "chip.no_dup_removal" : false,

    // Regular expression to filter out reads
    // Any read that matches with this reg-ex pattern will be removed from outputs
    "chip.regex_filter_reads" : "chrM",

    // Subsample reads (0: no subsampling)
    // Subsampled reads will be used for all downsteam analyses including peak-calling
    "chip.subsample_reads" : 0,
    "chip.ctl_subsample_reads" : 0,

    // Cross-correlation analysis
    // Subsample reads for cross-corr. analysis only (0: no subsampling)
    // Subsampled reads will be used for cross-corr. analysis only
    "chip.xcor_subsample_reads" : 15000000,

    // Keep irregular chromosome names
    // Use this for custom genomes without canonical chromosome names (chr1, chrX, ...)
    "chip.keep_irregular_chr_in_bfilt_peak" : false,

    // Choosing an appropriate control for each replicate
    // Always use a pooled control to compare with each replicate.
    // If a single control is given then use it.
    "chip.always_use_pooled_ctl" : false,
    // If ratio of depth between controls is higher than this
    // then always use a pooled control for all replicates.
    "chip.ctl_depth_ratio" : 1.2,

    // Cap number of peaks called from a peak-caller (MACS2)
    "chip.macs2_cap_num_peak" : 500000,
    // P-value threshold for MACS2 (macs2 callpeak -p)
    "chip.pval_thresh" : 0.01,

    // IDR (irreproducible discovery rate)
    // Threshold for IDR
    "chip.idr_thresh" : 0.05,

    // Cap number of peaks called from a peak-caller (SPP)
    "chip.spp_cap_num_peak" : 300000,

    ////////// 5) Resource settings //////////

    // Resources defined here are PER REPLICATE.
    // Therefore, total number of cores will be MAX(["chip.bwa_cpu"] x [NUMBER_OF_REPLICATES], "chip.spp_cpu" x 2 x [NUMBER_OF_REPLICATES])
    // because bwa and spp are bottlenecking tasks of the pipeline.
    // Use this total number of cores if you manually qsub or sbatch your job (using local mode of our pipeline).
    // "disks" is used for Google Cloud and DNANexus only.

    "chip.bwa_cpu" : 4,
    "chip.bwa_mem_mb" : 20000,
    "chip.bwa_time_hr" : 48,
    "chip.bwa_disks" : "local-disk 100 HDD",

    "chip.filter_cpu" : 2,
    "chip.filter_mem_mb" : 20000,
    "chip.filter_time_hr" : 24,
    "chip.filter_disks" : "local-disk 100 HDD",

    "chip.bam2ta_cpu" : 2,
    "chip.bam2ta_mem_mb" : 10000,
    "chip.bam2ta_time_hr" : 6,
    "chip.bam2ta_disks" : "local-disk 100 HDD",

    "chip.spr_mem_mb" : 16000,

    "chip.fingerprint_cpu" : 2,
    "chip.fingerprint_mem_mb" : 12000,
    "chip.fingerprint_time_hr" : 6,
    "chip.fingerprint_disks" : "local-disk 100 HDD",

    "chip.xcor_cpu" : 2,
    "chip.xcor_mem_mb" : 16000,
    "chip.xcor_time_hr" : 24,
    "chip.xcor_disks" : "local-disk 100 HDD",

    "chip.macs2_mem_mb" : 16000,
    "chip.macs2_time_hr" : 24,
    "chip.macs2_disks" : "local-disk 100 HDD",

    "chip.spp_cpu" : 2,
    "chip.spp_mem_mb" : 16000,
    "chip.spp_time_hr" : 72,
    "chip.spp_disks" : "local-disk 100 HDD",
}
```

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
|DNANexus (CLI, Azure)|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:data/pipeline-genome-data/[GENOME]_dx_azure.tsv`|
|DNANExus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data)|
|DNANExus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data)|
|Stanford Sherlock|`/home/groups/cherry/encode/pipeline_genome_data/[GENOME]_sherlock.tsv`|
|Stanford SCG|`/reference/ENCODE/pipeline_genome_data/[GENOME]_scg.tsv`|
|Local/SLURM/SGE|You need to [build a genome database](build_genome_database.md). |
