AQUAS Transcription Factor and Histone ChIP-Seq processing pipeline
===================================================

# Directories
* `backends/` : Backend configuration files (`.conf`)
* `workflow_opts/` : Workflow option files (`.json`)
* `examples/` : input JSON examples (SE and PE)
* `genome/` : genome data TSV files
* `src/` : Python script for each task in WDL
* `installers/` : dependency/genome data installers for systems (Local, SGE and SLURM) without docker support
* `docker_image/` : Dockerfile and MySQL DB initialization script

# Usage

See [Usage](https://github.com/encode-dcc/wdl-pipelines/blob/master/USAGE.md).

# Input JSON

Optional parameters and flags are marked with `?`. **`Input` in this document does not mean `control` here. `Control` is explicitly referred to as `control`.

1) Reference genome

    Currently supported genomes:

    * hg38: ENCODE [GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)
    * mm10: ENCODE [mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)
    * hg19: ENCODE [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)
    * mm9: [mm9, NCBI Build 37](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit)

    This TSV file has all genome specific data parameters and file path/URIs. Choose one of TSVs in `genome` directory.

    * `"chip.genome_tsv"` : TSV file path/URI.

2) Input genome data files
    Choose any genome data type you want to start with and do not define all others.
    
    * `"chip.fastqs"`? : 3-dimensional array with FASTQ file path/URI.
        - 1st dimension: replicate ID
        - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
        - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
    * `"chip.bams"`? : Array of raw (unfiltered) BAM file path/URI.
        - 1st dimension: replicate ID
    * `"chip.nodup_bams"`? : Array of filtered (deduped) BAM file path/URI.
        - 1st dimension: replicate ID
    * `"chip.tas"`? : Array of TAG-ALIGN file path/URI.
        - 1st dimension: replicate ID
    * `"chip.peaks"`? : Array of NARROWPEAK file path/URI.
        - 1st dimension: replicate ID
    * `"chip.peaks_pr1"`? : Array of NARROWPEAK file path/URI for 1st self pseudo replicate of replicate ID.
        - 1st dimension: replicate ID
    * `"chip.peaks_pr2"`? : Array of NARROWPEAK file path/URI for 2nd self pseudo replicate of replicate ID.
        - 1st dimension: replicate ID
    * `"chip.peak_ppr1"`? : NARROWPEAK file path/URI for pooled 1st pseudo replicates.
    * `"chip.peak_ppr2"`? : NARROWPEAK file path/URI for pooled 2nd pseudo replicates.
    * `"chip.peak_pooled"`? : NARROWPEAK file path/URI for pooled replicate.

    If starting from peaks then always define `"chip.peaks"`. Define `"chip.peaks_pr1"`, `"chip.peaks_pr2"`, `"chip.peak_pooled"`, `"chip.peak_ppr1"` and `"chip.peak_ppr2"` according to the following rules:    

    ```
    if num_rep>1:
        if true_rep_only: peak_pooled, 
        else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
    else:
        if true_rep_only: "not the case!"
        else: peaks_pr1[], peaks_pr2[]
    ```

    Default peak caller (`"chip.peak_caller"`) for TF (`"chip.pipeline_type":"tf"`) ChIP-Seq pipeline and Histone ChIP-Seq pipeline (`"chip.pipeline_type":"histone"`) are 'spp' and 'macs2', respectively. However you can also manually specify a peak caller for these pipeline types. `macs2` can work without controls but `spp` cannot. Therefore, if a peak caller is chosen as `spp` by default or by a workflow parameter then make sure to define the following control data files. Choose any genome data type you want to start with and do not define all others.

    * `"chip.ctl_fastqs"`? : 3-dimensional array with control FASTQ file path/URI.
        - 1st dimension: replicate ID
        - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
        - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
    * `"chip.ctl_bams"`? : Array of raw (unfiltered) control BAM file path/URI.
        - 1st dimension: replicate ID
    * `"chip.ctl_nodup_bams"`? : Array of filtered (deduped) control BAM file path/URI.
        - 1st dimension: replicate ID
    * `"chip.ctl_tas"`? : Array of control TAG-ALIGN file path/URI.
        - 1st dimension: replicate ID

    You can mix up different data types for IP replicates and controls. For example,
    ```
    input.json
    {
        "chip.paired_end" : false,
        "chip.bams" : ["rep1.bam","rep2.bam"],
        ...
        "chip.ctl_tas" : ["ctl1.tagAlign.gz","ctl2.tagAlign.gz"],
        ...
    }
    ```

3) Pipeline settings

    Pipeline type (chip-seq or DNase-Seq) : 'tf' ChIP-Seq always requires controls.

    * `"chip.pipeline_type` : `tf` for Transcription Factor ChIP-Seq. `histone` for Histone Factor ChIP-Seq.

    Input data endedness.

    * `"chip.paired_end"` : Set it as `true` if input dataset is paired end, otherwise `false`.

    Other important settings.

    * `"chip.align_only`? : Disable all downstream analysis after mapping.
    * `"chip.true_rep_only"`? : Set it as `true` to disable all analyses (including IDR, naive-overlap and reproducibility QC) related to pseudo replicates. This flag suppresses `"chip.enable_idr"`.

4) Trim FASTQ settings (**for paired end dataset only**)

    * `"chip.trim_fastq.trim_bp`? : For paired end dataset only. Number of basepairs after trimming FASTQ. It's 50 by default. Trimmed FASTQS is only used for cross-correlation analysis. FASTQ mapping is not affected by this parameter.

5) Filter/dedup (post-alignment) settings

    * `"chip.filter.dup_marker"`? : Dup marker. Choose between `picard` (default) and `sambamba`.
    * `"chip.filter.mapq_thresh"`? : Threshold for low MAPQ reads removal.
    * `"chip.filter.no_dup_removal"`? : No dup reads removal when filtering BAM.

6) BAM-2-TAGALIGN settings

    Pipeline filters out chrM reads by default.

    * `"chip.bam2ta.regex_grep_v_ta"`? : Perl-style regular expression pattern to remove matching reads from TAGALIGN (default: `chrM`).
    * `"chip.bam2ta.subsample"`? : Number of reads to subsample TAGALIGN. Subsampled TAGALIGN will be used for all downstream analysis (MACS2, IDR, naive-overlap).

7) Choose control settings

    * `"chip.choose_ctl.ctl_depth_ratio"`? : if ratio between controls is higher than this then always use pooled control for all exp rep.
    * `"chip.choose_ctl.always_use_pooled_ctl"`? : Always use pooled control for all exp replicates (ignoring `ctl_depth_ratio`).

8) Cross correlation analysis settings

    * `"chip.xcor.subsample"`? : Number of reads to subsample TAGALIGN. Only one end (R1) will be used for cross correlation analysis. This will not affect downstream analysis.

9) MACS2 settings

    **DO NOT DEFINE MACS2 PARAMETERS IN `"chip.macs2"` SCOPE**. All MACS2 parameters must be defined in `"chipseq"` scope.

    * `"chip.macs2_cap_num_peak"`? : Cap number of raw peaks called from MACS2.
    * `"chip.pval_thresh"`? : P-value threshold.

10) SPP settings

    **DO NOT DEFINE SPP PARAMETERS IN `"chip.spp"` SCOPE**. All SPP parameters must be defined in `"chipseq"` scope.

    * `"chip.spp_cap_num_peak"`? : Cap number of raw peaks called from MACS2.

11) IDR settings

    **DO NOT DEFINE IDR PARAMETERS IN `"chip.idr"` SCOPE**. All IDR parameters must be defined in `"chipseq"` scope.

    * `"chip.enable_idr"`? : Set it as `true` to enable IDR on raw peaks.
    * `"chip.idr_thresh"`? : IDR threshold.

12) Resources

    **RESOURCES DEFINED IN `input.json` ARE PER TASK**. For example, if you have FASTQs for 2 replicates (2 tasks) and set `cpu` for `bwa` task as 4 then total number of cpu cores to map FASTQs is 2 x 4 = 8.

    CPU (`cpu`), memory (`mem_mb`) settings are used for submitting jobs to cluster engines (SGE and SLURM) and Cloud platforms (Google Cloud Platform, AWS, ...). VM instance type on cloud platforms will be automatically chosen according to each task's `cpu` and `mem_mb`. Number of cores for tasks without `cpu` parameter is fixed at 1.

    * `"chip.merge_fastq.cpu"`? : Number of cores for `merge_fastq` (default: 2).
    * `"chip.bwa.cpu"`? : Number of cores for `bwa` (default: 4).
    * `"chip.filter.cpu"`? : Number of cores for `filter` (default: 2).
    * `"chip.bam2ta.cpu"`? : Number of cores for `bam2ta` (default: 2).
    * `"chip.xcor.cpu"`? : Number of cores for `xcor` (default: 2).
    * `"chip.spp_cpu"`? : Number of cores for `spp` (default: 2).
    * `"chip.merge_fastq.mem_mb"`? : Max. memory limit in MB for `merge_fastq` (default: 10000).
    * `"chip.bwa.mem_mb"`? : Max. memory limit in MB for `bwa` (default: 20000).
    * `"chip.filter.mem_mb"`? : Max. memory limit in MB for `filter` (default: 20000).
    * `"chip.bam2ta.mem_mb"`? : Max. memory limit in MB for `bam2ta` (default: 10000).
    * `"chip.spr.mem_mb"`? : Max. memory limit in MB for `spr` (default: 12000).
    * `"chip.xcor.mem_mb"`? : Max. memory limit in MB for `xcor` (default: 10000).
    * `"chip.macs2_mem_mb"`? : Max. memory limit in MB for `macs2` (default: 16000).
    * `"chip.spp_mem_mb"`? : Max. memory limit in MB for `spp` (default: 16000).

    Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

    * `"chip.merge_fastq.disks"`? : Disks for `merge_fastq` (default: "local-disk 100 HDD").
    * `"chip.bwa.disks"`? : Disks for `bwa` (default: "local-disk 100 HDD").
    * `"chip.filter.disks"`? : Disks for `filter` (default: "local-disk 100 HDD").
    * `"chip.bam2ta.disks"`? : Disks for `bam2ta` (default: "local-disk 100 HDD").
    * `"chip.xcor.disks"`? : Disks for `xcor` (default: "local-disk 100 HDD").
    * `"chip.spp_disks"`? : Disks for `spp` (default: "local-disk 100 HDD").
    * `"chip.macs2_disks"`? : Disks for `macs2` (default: "local-disk 100 HDD").

    Walltime (`time`) settings (for SGE and SLURM only).

    * `"chip.merge_fastq.time_hr"`? : Walltime for `merge_fastq` (default: 6).
    * `"chip.bwa.time_hr"`? : Walltime for `bwa` (default: 48).
    * `"chip.filter.time_hr"`? : Walltime for `filter` (default: 24).
    * `"chip.bam2ta.time_hr"`? : Walltime for `bam2ta` (default: 6).
    * `"chip.xcor.time_hr"`? : Walltime for `xcor` (default: 6).
    * `"chip.macs2_time_hr"`? : Walltime for `macs2` (default: 24).
    * `"chip.spp_time_hr"`? : Walltime for `spp` (default: 72).

13) QC report HTML/JSON

    * `"chip.qc_report.name"`? : Name of sample.
    * `"chip.qc_report.desc"`? : Description for sample.
