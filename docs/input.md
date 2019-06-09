# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file. We provide a set of template JSON files: [minimum](../examples/template.json) and [full](../examples/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

Please read through the following step-by-step instruction to compose a input JSON file.

## Pipeline metadata

Parameter|Description
---------|-----------
`chip.title`| Title for experiment which will be shown in a final HTML report
`chip.description`| Description for experiment which will be shown in a final HTML report

## Pipeline parameters

Parameter|Default|Description
---------|-------|-----------
`chip.pipeline_type`| `tf` | `tf` for TF ChIP-seq or `histone` for Histone ChIP-seq.
`chip.align_only`| false | Peak calling and its downstream analyses will be disabled. Useful if you just want to map your FASTQs into filtered BAMs/TAG-ALIGNs and don't want to call peaks on them.
`chip.true_rep_only` | false | Disable pseudo replicate generation and all related analyses

## Reference genome

All reference genome specific reference files/parameters can be defined in a single TSV file `chip.genome_tsv`. However, you can also individally define each file/parameter instead of a TSV file. If both a TSV file and individual parameters are defined, then individual parameters will override those defined in a TSV file. For example, if you define both `chip.genome_tsv` and `chip.blacklist`, then `chip.blacklist` will override that is defined in `chip.genome_tsv`. This is useful when you want to use your own for a specific parameter while keeping all the other parameters same as original.

Parameter|Type|Description
---------|-------|-----------
`chip.genome_tsv`| File | Choose one of the TSV files listed below or build your own
`chip.ref_fa`| File | Reference FASTA file
`chip.bwa_idx_tar`| File | BWA index TAR file (uncompressed) built from FASTA file with `bwa index`
`chip.chrsz`| File | 2-col chromosome sizes file built from FASTA file with `faidx`
`chip.blacklist`| File | 3-col BED file. Peaks overlapping these regions will be filtered out
`chip.gensz`| String | MACS2's genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)

We currently provide TSV files for 4 genomes as shown in the below table. `GENOME` should be `hg38`, `mm10`, `hg19` or `mm9`. You can [download/build](build_genome_database.md) it on your local computer. You can also [build a genome database for your own genome](build_genome_database.md).

Platform|Path/URI
-|-
Google Cloud Platform|`gs://encode-pipeline-genome-data/[GENOME]_google.tsv`
Stanford Sherlock|`/home/groups/cherry/encode/pipeline_genome_data/[GENOME]_sherlock.tsv`
Stanford SCG|`/reference/ENCODE/pipeline_genome_data/[GENOME]_scg.tsv`
Local/SLURM/SGE|You need to [build](build_genome_database.md) or [download]() a genome database]. 
DNAnexus (CLI)|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/[GENOME]_dx.tsv`
DNAnexus (CLI, Azure)|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/[GENOME]_dx_azure.tsv`
DNAnexus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data)
DNAnexus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data)

Additional information about each genome:

|Genome|Source|built from|
|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|

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
    "chip.nodup_bams" : [null, "rep2.bam", null, "rep4.nodup.bam", null],
    "chip.tas" : [null, null, null, null, "rep5.tagAlign.gz"]
}
```

## Optional filtering parameters

Parameter|Default|Description
---------|-------|-----------
`chip.mapq_thresh` | 30 | Threshold for mapped reads quality (samtools view -q)
`chip.dup_marker` | `picard` | Choose a dup marker between `picard` and `sambamba`. `picard` is recommended, use `sambamba` only when picard fails.
`chip.no_dup_removal` | false | Skip dup removal in a BAM filtering stage.

## Optional subsampling parameters

Parameter|Default|Description
---------|-------|-----------
`chip.subsample_reads` | 0 | Subsample reads (0: no subsampling). Subsampled reads will be used for all downsteam analyses including peak-calling
`chip.ctl_subsample_reads` | 0 | Subsample control reads. 
`chip.xcor_subsample_reads` | 15000000 | Subsample reads for cross-corr. analysis only (0: no subsampling). Subsampled reads will be used for cross-corr. analysis only

## Optional cross-correlation analysis parameters

Parameter|Default|Description
---------|-------|-----------
`chip.xcor_pe_trim_bp` | 50 | Trim R1 of paired ended fastqs for cross-correlation analysis only. Trimmed fastqs will not be used for any other analyses
`chip.use_filt_pe_ta_for_xcor` | false | Use filtered PE BAM/TAG-ALIGN for cross-correlation analysis ignoring the above trimmed R1 fastq

## Optional control parameters

Parameter|Default|Description
---------|-------|-----------
`chip.always_use_pooled_ctl` | false | Choosing an appropriate control for each replicate. Always use a pooled control to compare with each replicate. If a single control is given then use it.
`chip.ctl_depth_ratio` | 1.2 | If ratio of depth between controls is higher than this. then always use a pooled control for all replicates.

## Optional peak-calling parameters

Parameter|Default|Description
---------|-------|-----------
`chip.peak_caller`| `spp` for `tf` type<br>`macs2` for `histone` type| `spp` or `macs2`. `spp` requires control<br>`macs2` can work without controls
`chip.macs2_cap_num_peak` | 500000 | Cap number of peaks called from a peak-caller (MACS2)
`chip.pval_thresh` | 0.01 | P-value threshold for MACS2 (macs2 callpeak -p)
`chip.idr_thresh` | 0.05 | Threshold for IDR (irreproducible discovery rate)
`chip.spp_cap_num_peak` | 300000 | Cap number of peaks called from a peak-caller (SPP)

## Optional pipeline flags

Parameter|Default|Description
---------|-------|-----------
`chip.disable_fingerprint` | false | Disable deeptools fingerprint (JS distance)
`chip.enable_count_signal_track` | false | Enable count signal track generation
`chip.keep_irregular_chr_in_bfilt_peak` | false | Keep irregular chromosome names. Use this for custom genomes without canonical chromosome names (chr1, chrX, ...)

## Other optional parameters

Parameter|Default|Description
---------|-------|-----------
`chip.mito_chr_name` | `chrM` | Name of mito chromosome. THIS IS NOT A REG-EX! you can define only one chromosome name for mito.
`chip.regex_filter_reads` | `chrM` | Regular expression to filter out reads with given chromosome name (1st column of BED/TAG-ALIGN). Any read with chr name that matches with this reg-ex pattern will be removed from outputs If your have changed the above parameter `chip.mito_chr_name` and still want to filter out mito reads then make sure that `chip.mito_chr_name` and `chip.regex_filter_reads` are the same

## Resource parameters

> **WARNING*: It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task. The following parameters are provided for users who want to run our pipeline with Caper's `local` on HPCs and 2).

Resources defined here are PER REPLICATE. Therefore, total number of cores will be MAX(`chip.bwa_cpu` x `NUMBER_OF_REPLICATES`, `chip.spp_cpu` x 2 x `NUMBER_OF_REPLICATES`) because bwa and spp are bottlenecking tasks of the pipeline. Use this total number of cores if you manually `qsub` or `sbatch` your job (using local mode of Caper). `disks` is used for Google Cloud and DNAnexus only.

Parameter|Default
---------|-------
`chip.bwa_cpu` | 4
`chip.bwa_mem_mb` | 20000
`chip.bwa_time_hr` | 48
`chip.bwa_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`chip.filter_cpu` | 2
`chip.filter_mem_mb` | 20000
`chip.filter_time_hr` | 24
`chip.filter_disks` | `local-disk 100 HDD`

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
`chip.fingerprint_cpu` | 2
`chip.fingerprint_mem_mb` | 12000
`chip.fingerprint_time_hr` | 6
`chip.fingerprint_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`chip.xcor_cpu` | 2
`chip.xcor_mem_mb` | 16000
`chip.xcor_time_hr` | 24
`chip.xcor_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`chip.macs2_mem_mb` | 16000
`chip.macs2_time_hr` | 24
`chip.macs2_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`chip.spp_cpu` | 2
`chip.spp_mem_mb` | 16000
`chip.spp_time_hr` | 72
`chip.spp_disks` | `local-disk 100 HDD`
```
