version 1.0

workflow chip {
    String pipeline_ver = 'dev-v1.7.2'

    meta {
        version: 'dev-v1.7.2'
        author: 'Jin wook Lee (leepc12@gmail.com) at ENCODE-DCC'
        description: 'ENCODE TF/Histone ChIP-Seq pipeline'
        specification_document: 'https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit?usp=sharing'

        caper_docker: 'encodedcc/chip-seq-pipeline:dev-v1.7.2'
        caper_singularity: 'docker://encodedcc/chip-seq-pipeline:dev-v1.7.2'
        croo_out_def: 'https://storage.googleapis.com/encode-pipeline-output-definition/chip.croo.v5.json'

        parameter_group: {
            pipeline_metadata: {
                title: 'Pipeline metadata',
                description: 'Metadata for a pipeline (e.g. title and description).'
            },
            reference_genome: {
                title: 'Reference genome',
                description: 'Genome specific files. e.g. reference FASTA, bowtie2 index, chromosome sizes file.',
                help: 'Choose one chip.genome_tsv file that defines all genome specific parameters in it or define each genome specific parameter in input JSON to override those defined in genome TSV file. If you use Caper then use https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/[GENOME]_caper.tsv. Caper will automatically download/install all files defined in such TSV. Otherwise download genome TSV file by using a shell script (scripts/download_genome_data.sh [GENOME] [DEST_DIR]). Supported genomes are hg38, hg19, mm10 and mm9. See pipeline documentation if you want to build genome database from your own FASTA file. If some genome data are missing then analyses using such data will be skipped.'
            },
            input_genomic_data: {
                title: 'Input genomic data',
                description: 'Genomic input files for experiment.',
                help: 'Pipeline can start with any types of experiment data (e.g. FASTQ, BAM, NODUP_BAM, TAG-ALIGN, PEAK). Choose one type and leave others empty. FASTQs have a variable for each biological replicate. e.g. chip.fastqs_rep1_R1 and chip.fastqs_rep2_R1. You can define up to 10 experiment replicates. For other types, there is an array to define file for each biological replicate. e.g. chip.bams: ["rep1.bam", "rep1.bam"]. Define sequential endedness with chip.paired_end, if you have mixed SE and PE replicates then define chip.paired_ends instead for each replicate. e.g. chip.paired_ends: [false, true].'
            },
            input_genomic_data_control: {
                title: 'Input genomic data (control)',
                description: 'Genomic input files for control. TF ChIP-seq requires control for peak calling but histone ChIP-seq does not.',
                help: 'Pipeline can start with any types of control data (e.g. FASTQ, BAM, NODUP_BAM, TAG-ALIGN). Choose one type and leave others empty. FASTQs have a variable for each control replicate. e.g. chip.ctl_fastqs_rep1_R1 and chip.ctl_fastqs_rep2_R1. You can define up to 10 control replicates. For other types, there is an array to define file for each control replicate. e.g. chip.ctl_bams: ["ctl1.bam", "ctl1.bam"]. Define sequential endedness with chip.ctl_paired_end, if you have mixed SE and PE control replicates then define chip.ctl_paired_ends instead for each replicate. e.g. chip.ctl_paired_ends: [false, true]. If none of these are defined, pipeline will use chip.paired_end for controls.'
            },
            pipeline_parameter: {
                title: 'Pipeline parameter',
                description: 'Pipeline type and flags to turn on/off analyses.',
                help: 'Use chip.align_only to align FASTQs without peak calling.'
            },
            alignment: {
                title: 'Alignment',
                description: 'Parameters for alignment.',
                help: 'Pipeline can crop FASTQs (chip.crop_length > 0) with tolerance (chip.crop_length_tol) before mapping.'
            },
            peak_calling: {
                title: 'Peak calling',
                description: 'Parameters for peak calling.',
                help: 'This group includes statistical thresholds for peak-calling or post-peak-calling analyses: p-val, FDR, IDR. It also include parameters for control choosing/subsampling. All control replicates are pooled and pooled control is used for peak calling against each experiment replicate by default (see chip.always_use_pooled_ctl). Pipeline compares read depth of experiment replicate and a chosen control. It also compare read depth of controls. If control is too deep then it is subsampled.'
            },
            resource_parameter: {
                title: 'Resource parameter',
                description: 'Number of CPUs (threads), max. memory and walltime for tasks.',
                help: 'Resource settings are used for determining an instance type on cloud backends (e.g. GCP, AWS) and used for submitting tasks to a cluster engine (e.g. SLURM, SGE, ...). Walltime (chip.*_time_hr) is only used for cluster engines. Other tasks default to use 1 CPU and 4GB of memory.'
            }
        }
    }
    input {
        # group: pipeline_metadata
        String title = 'Untitled'
        String description = 'No description'

        # group: reference_genome
        File? genome_tsv
        String? genome_name
        File? ref_fa
        File? bwa_idx_tar
        File? bowtie2_idx_tar
        File? chrsz
        File? blacklist
        File? blacklist2
        String? mito_chr_name
        String? regex_bfilt_peak_chr_name
        String? gensz
        File? custom_aligner_idx_tar

        # group: input_genomic_data
        Boolean? paired_end
        Array[Boolean] paired_ends = []
        Array[File] fastqs_rep1_R1 = []
        Array[File] fastqs_rep1_R2 = []
        Array[File] fastqs_rep2_R1 = []
        Array[File] fastqs_rep2_R2 = []
        Array[File] fastqs_rep3_R1 = []
        Array[File] fastqs_rep3_R2 = []
        Array[File] fastqs_rep4_R1 = []
        Array[File] fastqs_rep4_R2 = []
        Array[File] fastqs_rep5_R1 = []
        Array[File] fastqs_rep5_R2 = []
        Array[File] fastqs_rep6_R1 = []
        Array[File] fastqs_rep6_R2 = []
        Array[File] fastqs_rep7_R1 = []
        Array[File] fastqs_rep7_R2 = []
        Array[File] fastqs_rep8_R1 = []
        Array[File] fastqs_rep8_R2 = []
        Array[File] fastqs_rep9_R1 = []
        Array[File] fastqs_rep9_R2 = []
        Array[File] fastqs_rep10_R1 = []
        Array[File] fastqs_rep10_R2 = []
        Array[File?] bams = []
        Array[File?] nodup_bams = []
        Array[File?] tas = []
        Array[File?] peaks = []
        Array[File?] peaks_pr1 = []
        Array[File?] peaks_pr2 = []
        File? peak_ppr1
        File? peak_ppr2
        File? peak_pooled

        Boolean? ctl_paired_end
        Array[Boolean] ctl_paired_ends = []
        Array[File] ctl_fastqs_rep1_R1 = []
        Array[File] ctl_fastqs_rep1_R2 = []
        Array[File] ctl_fastqs_rep2_R1 = []
        Array[File] ctl_fastqs_rep2_R2 = []
        Array[File] ctl_fastqs_rep3_R1 = []
        Array[File] ctl_fastqs_rep3_R2 = []
        Array[File] ctl_fastqs_rep4_R1 = []
        Array[File] ctl_fastqs_rep4_R2 = []
        Array[File] ctl_fastqs_rep5_R1 = []
        Array[File] ctl_fastqs_rep5_R2 = []
        Array[File] ctl_fastqs_rep6_R1 = []
        Array[File] ctl_fastqs_rep6_R2 = []
        Array[File] ctl_fastqs_rep7_R1 = []
        Array[File] ctl_fastqs_rep7_R2 = []
        Array[File] ctl_fastqs_rep8_R1 = []
        Array[File] ctl_fastqs_rep8_R2 = []
        Array[File] ctl_fastqs_rep9_R1 = []
        Array[File] ctl_fastqs_rep9_R2 = []
        Array[File] ctl_fastqs_rep10_R1 = []
        Array[File] ctl_fastqs_rep10_R2 = []
        Array[File?] ctl_bams = []
        Array[File?] ctl_nodup_bams = []
        Array[File?] ctl_tas = []

        # group: pipeline_parameter
        String pipeline_type
        Boolean align_only = false
        Boolean redact_nodup_bam = false
        Boolean true_rep_only = false
        Boolean enable_count_signal_track = false
        Boolean enable_jsd = true
        Boolean enable_gc_bias = true

        # group: alignment
        String aligner = 'bowtie2'
        File? custom_align_py
        Boolean use_bwa_mem_for_pe = false
        Int bwa_mem_read_len_limit = 70
        Boolean use_bowtie2_local_mode = false
        Int crop_length = 0
        Int crop_length_tol = 2
        String trimmomatic_phred_score_format = 'auto'
        Int xcor_trim_bp = 50
        Boolean use_filt_pe_ta_for_xcor = false
        String dup_marker = 'picard'
        Boolean no_dup_removal = false
        Int mapq_thresh = 30
        Array[String] filter_chrs = []
        Int subsample_reads = 0
        Int ctl_subsample_reads = 0
        Int xcor_subsample_reads = 15000000
        Int xcor_exclusion_range_min = -500
        Int? xcor_exclusion_range_max

        # group: peak_calling
        Int ctl_depth_limit = 200000000
        Float exp_ctl_depth_ratio_limit = 5.0
        Array[Int?] fraglen = []
        String? peak_caller
        Boolean always_use_pooled_ctl = true
        Float ctl_depth_ratio = 1.2
        Int? cap_num_peak
        Float pval_thresh = 0.01
        Float fdr_thresh = 0.01
        Float idr_thresh = 0.05

        # group: resource_parameter
        Int align_cpu = 6
        Float align_bowtie2_mem_factor = 0.15
        Float align_bwa_mem_factor = 1.0
        Int align_time_hr = 48
        Float align_bowtie2_disk_factor = 8.0
        Float align_bwa_disk_factor = 8.0

        Int filter_cpu = 4
        Float filter_mem_factor = 0.4
        Int filter_time_hr = 24
        Float filter_disk_factor = 8.0

        Int bam2ta_cpu = 2
        Float bam2ta_mem_factor = 0.35
        Int bam2ta_time_hr = 6
        Float bam2ta_disk_factor = 4.0

        Float spr_mem_factor = 13.5
        Float spr_disk_factor = 18.0

        Int jsd_cpu = 4
        Float jsd_mem_factor = 0.1
        Int jsd_time_hr = 6
        Float jsd_disk_factor = 2.0

        Int xcor_cpu = 2
        Float xcor_mem_factor = 1.0
        Int xcor_time_hr = 24
        Float xcor_disk_factor = 4.5

        Float subsample_ctl_mem_factor = 14.0
        Float subsample_ctl_disk_factor = 15.0

        Float macs2_signal_track_mem_factor = 12.0
        Int macs2_signal_track_time_hr = 24
        Float macs2_signal_track_disk_factor = 80.0

        Int call_peak_cpu = 6
        Float call_peak_spp_mem_factor = 5.0
        Float call_peak_macs2_mem_factor = 5.0
        Int call_peak_time_hr = 72
        Float call_peak_spp_disk_factor = 5.0
        Float call_peak_macs2_disk_factor = 30.0

        String? align_trimmomatic_java_heap
        String? filter_picard_java_heap
        String? gc_bias_picard_java_heap
    }

    parameter_meta {
        title: {
            description: 'Experiment title.',
            group: 'pipeline_metadata',
            example: 'ENCSR936XTK (subsampled 1/50)'
        }
        description: {
            description: 'Experiment description.',
            group: 'pipeline_metadata',
            example: 'ZNF143 ChIP-seq on human GM12878 (subsampled 1/50)'
        }
        genome_tsv: {
            description: 'Reference genome database TSV.',
            group: 'reference_genome',
            help: 'This TSV files includes all genome specific parameters (e.g. reference FASTA, bowtie2 index). You can still invidiaully define any parameters in it. Parameters defined in input JSON will override those defined in genome TSV.',
            example: 'https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/hg38_caper.tsv'
        }
        genome_name: {
            description: 'Genome name.',
            group: 'reference_genome'
        }
        ref_fa: {
            description: 'Reference FASTA file.',
            group: 'reference_genome'
        }
        bowtie2_idx_tar: {
            description: 'BWA index TAR file.',
            group: 'reference_genome'
        }
        custom_aligner_idx_tar: {
            description: 'Index TAR file for a custom aligner. To use a custom aligner, define "chip.custom_align_py" too.',
            group: 'reference_genome'
        }
        chrsz: {
            description: '2-col chromosome sizes file.',
            group: 'reference_genome'
        }
        blacklist: {
            description: 'Blacklist file in BED format.',
            group: 'reference_genome',
            help: 'Peaks will be filtered with this file.'
        }
        blacklist2: {
            description: 'Secondary blacklist file in BED format.',
            group: 'reference_genome',
            help: 'If it is defined, it will be merged with chip.blacklist. Peaks will be filtered with merged blacklist.'
        }
        mito_chr_name: {
            description: 'Mitochondrial chromosome name.',
            group: 'reference_genome',
            help: 'e.g. chrM, MT. Mitochondrial reads defined here will be filtered out during filtering BAMs in "filter" task.'
        }
        regex_bfilt_peak_chr_name: {
            description: 'Reg-ex for chromosomes to keep while filtering peaks.',
            group: 'reference_genome',
            help: 'Chromosomes defined here will be kept. All other chromosomes will be filtered out in .bfilt. peak file. This is done along with blacklist filtering peak file.'
        }
        gensz: {
            description: 'Genome sizes. "hs" for human, "mm" for mouse or sum of 2nd columnin chromosome sizes file.',
            group: 'reference_genome'
        }
        paired_end: {
            description: 'Sequencing endedness.',
            group: 'input_genomic_data',
            help: 'Setting this on means that all replicates are paired ended. For mixed samples, use chip.paired_ends array instead.',
            example: true
        }
        paired_ends: {
            description: 'Sequencing endedness array (for mixed SE/PE datasets).',
            group: 'input_genomic_data',
            help: 'Whether each biological replicate is paired ended or not.'
        }
        fastqs_rep1_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from FASTQs files. Pipeline can start from any type of inputs (e.g. FASTQs, BAMs, ...). Choose one type and fill paramters for that type and leave other undefined. Especially for FASTQs, we have individual variable for each biological replicate to allow FASTQs of technical replicates can be merged. Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep1_R2). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/rep1-R1.subsampled.50.fastq.gz'
            ]
        }
        fastqs_rep1_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep1_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/rep1-R2.subsampled.50.fastq.gz'
            ]
        }
        fastqs_rep2_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 2.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep2_R2). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/rep2-R1.subsampled.50.fastq.gz'
            ]
        }
        fastqs_rep2_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 2.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep2_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/rep2-R2.subsampled.50.fastq.gz'
            ]
        }
        fastqs_rep3_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 3.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep3_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep3_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 3.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep3_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep4_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 4.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep4_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep4_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 4.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep4_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep5_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 5.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep5_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep5_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 5.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep5_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep6_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 6.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep6_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep6_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 6.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep6_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep7_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 7.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep7_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep7_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 7.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep7_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep8_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 8.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep8_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep8_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 8.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep8_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep9_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 9.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep9_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep9_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 9.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep9_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep10_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 10.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.fastqs_rep10_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep10_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 10.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.fastqs_rep10_R1). These FASTQs are usually technical replicates to be merged.'
        }
        bams: {
            description: 'List of unfiltered/raw BAM files for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from BAM files. Unfiltered/raw BAM file generated from aligner (e.g. bowtie2). Each entry for each biological replicate. e.g. [rep1.bam, rep2.bam, rep3.bam, ...].'
        }
        nodup_bams: {
            description: 'List of filtered/deduped BAM files for each biological replicate',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from filtered BAM files. Filtered/deduped BAM file. Each entry for each biological replicate. e.g. [rep1.nodup.bam, rep2.nodup.bam, rep3.nodup.bam, ...].'
        }
        tas: {
            description: 'List of TAG-ALIGN files for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from TAG-ALIGN files. TAG-ALIGN is in a 6-col BED format. It is a simplified version of BAM. Each entry for each biological replicate. e.g. [rep1.tagAlign.gz, rep2.tagAlign.gz, ...].'
        }
        peaks: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Each entry for each biological replicate. e.g. [rep1.narrowPeak.gz, rep2.narrowPeak.gz, ...]. Define other PEAK parameters (e.g. chip.peaks_pr1, chip.peak_pooled) according to your flag settings (e.g. chip.true_rep_only) and number of replicates. If you have more than one replicate then define chip.peak_pooled, chip.peak_ppr1 and chip.peak_ppr2. If chip.true_rep_only flag is on then do not define any parameters (chip.peaks_pr1, chip.peaks_pr2, chip.peak_ppr1 and chip.peak_ppr2) related to pseudo replicates.'
        }
        peaks_pr1: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for pseudo-replicate 1 of each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if chip.true_rep_only flag is off.'
        }
        peaks_pr2: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for pseudo-replicate 2 of each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if chip.true_rep_only flag is off.'
        }
        peak_pooled: {
            description: 'NARROWPEAK file for pooled true replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if you have multiple biological replicates. Pooled true replicate means analysis on pooled biological replicates.'
        }
        peak_ppr1: {
            description: 'NARROWPEAK file for pooled pseudo replicate 1.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if you have multiple biological replicates and chip.true_rep_only flag is off. PPR1 means analysis on pooled 1st pseudo replicates. Each biological replicate is shuf/split into two pseudos. This is a pooling of each replicate\'s 1st pseudos.'
        }
        peak_ppr2: {
            description: 'NARROWPEAK file for pooled pseudo replicate 2.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if you have multiple biological replicates and chip.true_rep_only flag is off. PPR1 means analysis on pooled 2nd pseudo replicates. Each biological replicate is shuf/split into two pseudos. This is a pooling of each replicate\'s 2nd pseudos.'
        }

        ctl_paired_end: {
            description: 'Sequencing endedness for all controls.',
            group: 'input_genomic_data_control',
            help: 'Setting this on means that all control replicates are paired ended. For mixed controls, use chip.ctl_paired_ends array instead.'
        }
        ctl_paired_ends: {
            description: 'Sequencing endedness array for mixed SE/PE controls.',
            group: 'input_genomic_data_control',
            help: 'Whether each control replicate is paired ended or not.'
        }
        ctl_fastqs_rep1_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 1.',
            group: 'input_genomic_data_control',
            help: 'Define if you want to start pipeline from FASTQs files. Pipeline can start from any type of controls (e.g. FASTQs, BAMs, ...). Choose one type and fill paramters for that type and leave other undefined.  Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep1_R2).',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/ctl1-R1.subsampled.80.fastq.gz'
            ]
        }
        ctl_fastqs_rep1_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 1.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep1_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/ctl1-R2.subsampled.80.fastq.gz'
            ]
        }
        ctl_fastqs_rep2_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 2.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep2_R2). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/ctl2-R1.subsampled.80.fastq.gz'
            ]
        }
        ctl_fastqs_rep2_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 2.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep2_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                'https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/fastq_subsampled/ctl2-R2.subsampled.80.fastq.gz'
            ]
        }
        ctl_fastqs_rep3_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 3.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep3_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep3_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 3.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep3_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep4_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 4.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep4_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep4_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 4.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep4_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep5_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 5.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep5_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep5_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 5.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep5_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep6_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 6.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep6_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep6_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 6.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep6_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep7_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 7.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep7_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep7_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 7.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep7_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep8_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 8.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep8_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep8_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 8.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep8_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep9_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 9.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep9_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep9_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 9.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep9_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep10_R1: {
            description: 'Read1 FASTQs to be merged for a control replicate 10.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read2 FASTQs (chip.ctl_fastqs_rep10_R2). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_fastqs_rep10_R2: {
            description: 'Read2 FASTQs to be merged for a control replicate 10.',
            group: 'input_genomic_data_control',
            help: 'Make sure that they are consistent with read1 FASTQs (chip.ctl_fastqs_rep10_R1). These FASTQs are usually technical replicates to be merged.'
        }
        ctl_bams: {
            description: 'List of unfiltered/raw BAM files for each control replicate.',
            group: 'input_genomic_data_control',
            help: 'Define if you want to start pipeline from BAM files. Unfiltered/raw BAM file generated from aligner (e.g. bowtie2). Each entry for each control replicate. e.g. [ctl1.bam, ctl2.bam, ctl3.bam, ...].'
        }
        ctl_nodup_bams: {
            description: 'List of filtered/deduped BAM files for each control replicate',
            group: 'input_genomic_data_control',
            help: 'Define if you want to start pipeline from filtered BAM files. Filtered/deduped BAM file. Each entry for each control replicate. e.g. [ctl1.nodup.bam, ctl2.nodup.bam, ctl3.nodup.bam, ...].'
        }
        ctl_tas: {
            description: 'List of TAG-ALIGN files for each biological replicate.',
            group: 'input_genomic_data_control',
            help: 'Define if you want to start pipeline from TAG-ALIGN files. TAG-ALIGN is in a 6-col BED format. It is a simplified version of BAM. Each entry for each control replicate. e.g. [ctl1.tagAlign.gz, ctl2.tagAlign.gz, ...].'
        }

        pipeline_type: {
            description: 'Pipeline type. tf for TF ChIP-Seq, histone for Histone ChIP-Seq or control for mapping controls only.',
            group: 'pipeline_parameter',
            help: 'Default peak caller is different for each type. spp For TF ChIP-Seq and macs2 for histone ChIP-Seq. Regardless of pipeline type, spp always requires controls but macs2 doesn\'t. For control mode, chip.align_only is automatically turned on and cross-correlation analysis is disabled. Do not define ctl_* for control mode. Define fastqs_repX_RY instead.',
            choices: ['tf', 'histone', 'control'],
            example: 'tf'
        }
        redact_nodup_bam: {
            description: 'Redact filtered/nodup BAM.',
            group: 'pipeline_parameter',
            help: 'Redact filtered/nodup BAM at the end of the filtering step (task filter). Raw BAM from the aligner (task align) will still remain unredacted. Quality metrics on filtered BAM will be calculated before being redacted. However, all downstream analyses (e.g. peak-calling) will be done on the redacted BAM. If you start from nodup BAM then this flag will not be active.'
        }
        align_only: {
            description: 'Align only mode.',
            group: 'pipeline_parameter',
            help: 'Reads will be aligned but there will be no peak-calling on them. It is turned on automatically if chip.pipeline_type is control.'
        }
        true_rep_only: {
            description: 'Disables all analyses related to pseudo-replicates.',
            group: 'pipeline_parameter',
            help: 'Pipeline generates 2 pseudo-replicate from one biological replicate. This flag turns off all analyses related to pseudos (with prefix/suffix pr, ppr).'
        }
        enable_count_signal_track: {
            description: 'Enables generation of count signal tracks.',
            group: 'pipeline_parameter'
        }
        enable_jsd: {
            description: 'Enables Jensen-Shannon Distance (JSD) plot generation.',
            group: 'pipeline_parameter'
        }
        enable_gc_bias: {
            description: 'Enables GC bias calculation.',
            group: 'pipeline_parameter'
        }

        aligner: {
            description: 'Aligner. bowtie2, bwa or custom',
            group: 'alignment',
            help: 'It is bowtie2 by default. To use a custom aligner, define chip.custom_align_py and chip.custom_aligner_idx_tar.',
            choices: ['bowtie2', 'bwa', 'custom'],
            example: 'bowtie2'
        }
        custom_align_py: {
            description: 'Python script for a custom aligner.',
            group: 'alignment',
            help: 'There is a template included in the documentation for inputs. Defining this parameter will automatically change "chip.aligner" to "custom". You should also define "chip.custom_aligner_idx_tar".'
        }
        use_bwa_mem_for_pe: {
            description: 'For paired end dataset with read length >= chip.bwa_mem_read_len_limit (default 70) bp, use bwa mem instead of bwa aln.',
            group: 'alignment',
            help: 'Use it only for paired end reads >= chip.bwa_mem_read_len_limit (default 70) bp. Otherwise keep using bwa aln.'
        }
        bwa_mem_read_len_limit {
            description: 'Read length limit for bwa mem (for PE FASTQs only).',
            group: 'alignment',
            help: 'If chip.use_bwa_mem_for_pe is activated and reads are shorter than this limit, then bwa aln will be used instead of bwa mem.'
        }
        use_bowtie2_local_mode: {
            description: 'Use bowtie2\'s local mode (soft-clipping).',
            group: 'alignment',
            help: 'This will add --local to bowtie2 command line so that it will replace the default end-to-end mode.'
        }
        crop_length: {
            description: 'Crop FASTQs\' reads longer than this length.',
            group: 'alignment',
            help: 'Also drop all reads shorter than chip.crop_length - chip.crop_length_tol.'
        }
        crop_length_tol: {
            description: 'Tolerance for cropping reads in FASTQs.',
            group: 'alignment',
            help: 'Drop all reads shorter than chip.crop_length - chip.crop_length_tol. Activated only when chip.crop_length is defined.'
        }
        trimmomatic_phred_score_format: {
            description: 'Base encoding (format) for Phred score in FASTQs.',
            group: 'alignment',
            choices: ['auto', 'phred33', 'phred64'],
            help: 'This is used for Trimmomatic only. It is auto by default, which means that Trimmomatic automatically detect it from FASTQs. Otherwise -phred33 or -phred64 will be passed to the Trimmomatic command line. Use this if you see an error like "Error: Unable to detect quality encoding".'
        }
        xcor_trim_bp: {
            description: 'Trim experiment read1 FASTQ (for both SE and PE) for cross-correlation analysis.',
            group: 'alignment',
            help: 'This does not affect alignment of experimental/control replicates. Pipeline additionaly aligns R1 FASTQ only for cross-correlation analysis only. This parameter is used for it.'
        }
        use_filt_pe_ta_for_xcor: {
            description: 'Use filtered PE BAM for cross-correlation analysis.',
            group: 'alignment',
            help: 'If not defined, pipeline uses SE BAM generated from trimmed read1 FASTQ for cross-correlation analysis.'
        }
        dup_marker: {
            description: 'Marker for duplicate reads. picard or sambamba.',
            group: 'alignment',
            help: 'picard for Picard MarkDuplicates or sambamba for sambamba markdup.',
            choices: ['picard', 'sambamba'],
            example: 'picard'
        }
        no_dup_removal: {
            description: 'Disable removal of duplicate reads during filtering BAM.',
            group: 'alignment',
            help: 'Duplicate reads are filtererd out during filtering BAMs to gerenate NODUP_BAM. This flag will keep all duplicate reads in NODUP_BAM. This flag does not affect naming of NODUP_BAM. NODUP_BAM will still have .nodup. suffix in its filename.'
        }
        mapq_thresh: {
            description: 'Threshold for low MAPQ reads removal.',
            group: 'alignment',
            help: 'Low MAPQ reads are filtered out while filtering BAM.'
        }
        filter_chrs: {
            description: 'List of chromosomes to be filtered out while filtering BAM.',
            group: 'alignment',
            help: 'It is empty by default, hence no filtering out of specfic chromosomes. It is case-sensitive. Use exact word for chromosome names.'
        }
        subsample_reads: {
            description: 'Subsample reads. Shuffle and subsample reads.',
            group: 'alignment',
            help: 'This affects all downstream analyses after filtering experiment BAM. (e.g. all TAG-ALIGN files, peak-calling). Reads will be shuffled only if actual number of reads in BAM exceeds this number. 0 means disabled.'
        }
        ctl_subsample_reads: {
            description: 'Subsample control reads. Shuffle and subsample control reads.',
            group: 'alignment',
            help: 'This affects all downstream analyses after filtering control BAM. (e.g. all TAG-ALIGN files, peak-calling). Reads will be shuffled only if actual number of reads in BAM exceeds this number. 0 means disabled.'
        }
        xcor_subsample_reads: {
            description: 'Subsample reads for cross-corrlelation analysis only.',
            group: 'alignment',
            help: 'This does not affect downstream analyses after filtering BAM. It is for cross-correlation analysis only.  0 means disabled.'
        }
        xcor_exclusion_range_min: {
            description: 'Exclusion minimum for cross-correlation analysis.',
            group: 'alignment',
            help: 'For run_spp.R -s. Make sure that it is consistent with default strand shift -s=-500:5:1500 in run_spp.R.'
        }
        xcor_exclusion_range_max: {
            description: 'Exclusion maximum for cross-coorrelation analysis.',
            group: 'alignment',
            help: 'For run_spp.R -s. If not defined default value of `max(read length + 10, 50)` for TF and `max(read_len + 10, 100)` for histone are used'
        }

        ctl_depth_limit: {
            description: 'Hard limit for chosen control\'s depth.',
            group: 'peak_calling',
            help: 'If control chosen by chip.always_use_pooled_ctl and chip.ctl_depth_ratio is deeper than this hard limit, then such control is subsampled.'
        }
        exp_ctl_depth_ratio_limit: {
            description: 'Second limit for chosen control\'s depth.',
            group: 'peak_calling',
            help: 'If control chosen by chip.always_use_pooled_ctl and chip.ctl_depth_ratio is deeper than experiment replicate\'s read depth multiplied by this factor then such control is subsampled down to maximum of multiplied value and hard limit chip.ctl_depth_limit.'
        }
        fraglen: {
            description: 'Fragment length for each biological replicate.',
            group: 'peak_calling',
            help: 'Fragment length is estimated by cross-correlation analysis, which is valid only when pipeline started from FASTQs. If defined, fragment length estimated by cross-correlation analysis is ignored.'
        }
        peak_caller: {
            description: 'Peak caller.',
            group: 'peak_calling',
            help: 'It is spp and macs2 by default for TF ChIP-seq and histone ChIP-seq, respectively. e.g. you can use macs2 for TF ChIP-Seq even though spp is by default for TF ChIP-Seq (chip.pipeline_type == tf).',
            choices: ['spp', 'macs2'],
            example: 'spp'
        }
        always_use_pooled_ctl: {
            description: 'Always choose a pooled control for each experiment replicate.',
            group: 'peak_calling',
            help: 'If turned on, ignores chip.ctl_depth_ratio.'
        }
        ctl_depth_ratio: {
            description: 'Maximum depth ratio between control replicates.',
            group: 'peak_calling',
            help: 'If ratio of depth between any two controls is higher than this, then always use a pooled control for all experiment replicates.'
        }

        cap_num_peak: {
            description: 'Upper limit on the number of peaks.',
            group: 'peak_calling',
            help: 'It is 30000000 and 50000000 by default for spp and macs2, respectively.'
        }
        pval_thresh: {
            description: 'p-value Threshold for MACS2 peak caller.',
            group: 'peak_calling',
            help: 'macs2 callpeak -p'
        }
        fdr_thresh: {
            description: 'FDR threshold for spp peak caller (phantompeakqualtools).',
            group: 'peak_calling',
            help: 'run_spp.R -fdr='
        }
        idr_thresh: {
            description: 'IDR threshold.',
            group: 'peak_calling'
        }

        align_cpu: {
            description: 'Number of cores for task align.',
            group: 'resource_parameter',
            help: 'Task align merges/crops/maps FASTQs.'
        }
        align_bowtie2_mem_factor: {
            description: 'Multiplication factor to determine memory required for task align with bowtie2 (default) as aligner.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        align_bwa_mem_factor: {
            description: 'Multiplication factor to determine memory required for task align with bwa as aligner.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        align_time_hr: {
            description: 'Walltime (h) required for task align.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        align_bowtie2_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task align with bowtie2 (default) as aligner.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required disk size of instance on GCP/AWS.'
        }
        align_bwa_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task align with bwa as aligner.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required disk size of instance on GCP/AWS.'
        }
        filter_cpu: {
            description: 'Number of cores for task filter.',
            group: 'resource_parameter',
            help: 'Task filter filters raw/unfiltered BAM to get filtered/deduped BAM.'
        }
        filter_mem_factor: {
            description: 'Multiplication factor to determine memory required for task filter.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        filter_time_hr: {
            description: 'Walltime (h) required for task filter.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        filter_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task filter.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required disk size of instance on GCP/AWS.'
        }
        bam2ta_cpu: {
            description: 'Number of cores for task bam2ta.',
            group: 'resource_parameter',
            help: 'Task bam2ta converts filtered/deduped BAM in to TAG-ALIGN (6-col BED) format.'
        }
        bam2ta_mem_factor: {
            description: 'Multiplication factor to determine memory required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        bam2ta_time_hr: {
            description: 'Walltime (h) required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        bam2ta_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task bam2ta.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        spr_mem_factor: {
            description: 'Multiplication factor to determine memory required for task spr.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        spr_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task spr.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        jsd_cpu: {
            description: 'Number of cores for task jsd.',
            group: 'resource_parameter',
            help: 'Task jsd plots Jensen-Shannon distance and metrics related to it.'
        }
        jsd_mem_factor: {
            description: 'Multiplication factor to determine memory required for task jsd.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        jsd_time_hr: {
            description: 'Walltime (h) required for task jsd.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        jsd_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task jsd.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        xcor_cpu: {
            description: 'Number of cores for task xcor.',
            group: 'resource_parameter',
            help: 'Task xcor does cross-correlation analysis (including a plot) on subsampled TAG-ALIGNs.'
        }
        xcor_mem_factor: {
            description: 'Multiplication factor to determine memory required for task xcor.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        xcor_time_hr: {
            description: 'Walltime (h) required for task xcor.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        xcor_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task xcor.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        subsample_ctl_mem_factor: {
            description: 'Multiplication factor to determine memory required for task subsample_ctl.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        subsample_ctl_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task subsample_ctl.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        call_peak_cpu: {
            description: 'Number of cores for task call_peak. IF MACS2 is chosen as peak_caller (or chip.pipeline_type is histone), then cpu will be fixed at 2.',
            group: 'resource_parameter',
            help: 'Task call_peak call peaks on TAG-ALIGNs by using SPP/MACS2 peak caller. MACS2 is single-threaded so cpu will be fixed at 2 for MACS2.'
        }
        call_peak_spp_mem_factor: {
            description: 'Multiplication factor to determine memory required for task call_peak with spp as peak_caller.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        call_peak_macs2_mem_factor: {
            description: 'Multiplication factor to determine memory required for task call_peak with macs2 as peak_caller.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        call_peak_time_hr: {
            description: 'Walltime (h) required for task call_peak.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        call_peak_spp_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task call_peak with spp as peak_caller.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        call_peak_macs2_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task call_peak with macs2 as peak_caller.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        macs2_signal_track_mem_factor: {
            description: 'Multiplication factor to determine memory required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        macs2_signal_track_time_hr: {
            description: 'Walltime (h) required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        macs2_signal_track_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        align_trimmomatic_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task align.',
            group: 'resource_parameter',
            help: 'Maximum memory for Trimmomatic. If not defined, 90% of align task\'s memory will be used.'
        }
        filter_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task filter.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools MarkDuplicates. If not defined, 90% of filter task\'s memory will be used.'
        }
        gc_bias_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task gc_bias.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectGcBiasMetrics. If not defined, 90% of gc_bias task\'s memory will be used.'
        }
    }

    # read genome data and paths
    if ( defined(genome_tsv) ) {
        call read_genome_tsv { input: genome_tsv = genome_tsv }
    }
    File ref_fa_ = select_first([ref_fa, read_genome_tsv.ref_fa])
    File? bwa_idx_tar_ = if defined(bwa_idx_tar) then bwa_idx_tar
        else read_genome_tsv.bwa_idx_tar
    File bowtie2_idx_tar_ = select_first([bowtie2_idx_tar, read_genome_tsv.bowtie2_idx_tar])
    File chrsz_ = select_first([chrsz, read_genome_tsv.chrsz])
    String gensz_ = select_first([gensz, read_genome_tsv.gensz])
    File? blacklist1_ = if defined(blacklist) then blacklist
        else read_genome_tsv.blacklist
    File? blacklist2_ = if defined(blacklist2) then blacklist2
        else read_genome_tsv.blacklist2
    # merge multiple blacklists
    # two blacklists can have different number of columns (3 vs 6)
    # so we limit merged blacklist's columns to 3
    Array[File] blacklists = select_all([blacklist1_, blacklist2_])
    if ( length(blacklists) > 1 ) {
        call pool_ta as pool_blacklist { input:
            tas = blacklists,
            col = 3,
        }
    }
    File? blacklist_ = if length(blacklists) > 1 then pool_blacklist.ta_pooled
        else if length(blacklists) > 0 then blacklists[0]
        else blacklist2_
    String mito_chr_name_ = select_first([mito_chr_name, read_genome_tsv.mito_chr_name])
    String regex_bfilt_peak_chr_name_ = select_first([regex_bfilt_peak_chr_name, read_genome_tsv.regex_bfilt_peak_chr_name])
    String genome_name_ = select_first([genome_name, read_genome_tsv.genome_name, basename(chrsz_)])

    ### temp vars (do not define these)
    String aligner_ = if defined(custom_align_py) then 'custom' else aligner
    String peak_caller_ = if pipeline_type=='tf' then select_first([peak_caller, 'spp'])
                        else select_first([peak_caller, 'macs2'])
    String peak_type_ = if peak_caller_=='spp' then 'regionPeak'
                        else 'narrowPeak'
    Boolean enable_idr = pipeline_type=='tf' # enable_idr for TF chipseq only
    String idr_rank_ = if peak_caller_=='spp' then 'signal.value'
                        else if peak_caller_=='macs2' then 'p.value'
                        else 'p.value'
    Int cap_num_peak_spp = 300000
    Int cap_num_peak_macs2 = 500000                        
    Int cap_num_peak_ = if peak_caller_ == 'spp' then select_first([cap_num_peak, cap_num_peak_spp])
        else select_first([cap_num_peak, cap_num_peak_macs2])
    Int mapq_thresh_ = mapq_thresh
    Boolean enable_xcor_ = if pipeline_type=='control' then false else true
    Boolean enable_count_signal_track_ = if pipeline_type=='control' then false else enable_count_signal_track
    Boolean enable_jsd_ = if pipeline_type=='control' then false else enable_jsd
    Boolean enable_gc_bias_ = if pipeline_type=='control' then false else enable_gc_bias
    Boolean align_only_ = if pipeline_type=='control' then true else align_only

    Float align_mem_factor_ = if aligner_ =='bowtie2' then align_bowtie2_mem_factor
        else align_bwa_mem_factor
    Float align_disk_factor_ = if aligner_ =='bowtie2' then align_bowtie2_disk_factor
        else align_bwa_disk_factor
    Float call_peak_mem_factor_ = if peak_caller_ =='spp' then call_peak_spp_mem_factor
        else call_peak_macs2_mem_factor
    Float call_peak_disk_factor_ = if peak_caller_ =='spp' then call_peak_spp_disk_factor
        else call_peak_macs2_disk_factor

    # temporary 2-dim fastqs array [rep_id][merge_id]
    Array[Array[File]] fastqs_R1 = 
        if length(fastqs_rep10_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1, fastqs_rep9_R1, fastqs_rep10_R1]
        else if length(fastqs_rep9_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1, fastqs_rep9_R1]
        else if length(fastqs_rep8_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1]
        else if length(fastqs_rep7_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1]
        else if length(fastqs_rep6_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1]
        else if length(fastqs_rep5_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1]
        else if length(fastqs_rep4_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1]
        else if length(fastqs_rep3_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1]
        else if length(fastqs_rep2_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1]
        else if length(fastqs_rep1_R1)>0 then
            [fastqs_rep1_R1]
        else []
    # no need to do that for R2 (R1 array will be used to determine presense of fastq for each rep)
    Array[Array[File]] fastqs_R2 = 
        [fastqs_rep1_R2, fastqs_rep2_R2, fastqs_rep3_R2, fastqs_rep4_R2, fastqs_rep5_R2,
        fastqs_rep6_R2, fastqs_rep7_R2, fastqs_rep8_R2, fastqs_rep9_R2, fastqs_rep10_R2]

    # temporary 2-dim ctl fastqs array [rep_id][merge_id]
    Array[Array[File]] ctl_fastqs_R1 = 
        if length(ctl_fastqs_rep10_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1, ctl_fastqs_rep5_R1,
            ctl_fastqs_rep6_R1, ctl_fastqs_rep7_R1, ctl_fastqs_rep8_R1, ctl_fastqs_rep9_R1, ctl_fastqs_rep10_R1]
        else if length(ctl_fastqs_rep9_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1, ctl_fastqs_rep5_R1,
            ctl_fastqs_rep6_R1, ctl_fastqs_rep7_R1, ctl_fastqs_rep8_R1, ctl_fastqs_rep9_R1]
        else if length(ctl_fastqs_rep8_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1, ctl_fastqs_rep5_R1,
            ctl_fastqs_rep6_R1, ctl_fastqs_rep7_R1, ctl_fastqs_rep8_R1]
        else if length(ctl_fastqs_rep7_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1, ctl_fastqs_rep5_R1,
            ctl_fastqs_rep6_R1, ctl_fastqs_rep7_R1]
        else if length(ctl_fastqs_rep6_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1, ctl_fastqs_rep5_R1,
            ctl_fastqs_rep6_R1]
        else if length(ctl_fastqs_rep5_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1, ctl_fastqs_rep5_R1]
        else if length(ctl_fastqs_rep4_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1, ctl_fastqs_rep4_R1]
        else if length(ctl_fastqs_rep3_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1, ctl_fastqs_rep3_R1]
        else if length(ctl_fastqs_rep2_R1)>0 then
            [ctl_fastqs_rep1_R1, ctl_fastqs_rep2_R1]
        else if length(ctl_fastqs_rep1_R1)>0 then
            [ctl_fastqs_rep1_R1]
        else []
    # no need to do that for R2 (R1 array will be used to determine presense of fastq for each rep)
    Array[Array[File]] ctl_fastqs_R2 = 
        [ctl_fastqs_rep1_R2, ctl_fastqs_rep2_R2, ctl_fastqs_rep3_R2, ctl_fastqs_rep4_R2, ctl_fastqs_rep5_R2,
        ctl_fastqs_rep6_R2, ctl_fastqs_rep7_R2, ctl_fastqs_rep8_R2, ctl_fastqs_rep9_R2, ctl_fastqs_rep10_R2]

    # temporary variables to get number of replicates
    #     WDLic implementation of max(A,B,C,...)
    Int num_rep_fastq = length(fastqs_R1)
    Int num_rep_bam = if length(bams)<num_rep_fastq then num_rep_fastq
        else length(bams)
    Int num_rep_nodup_bam = if length(nodup_bams)<num_rep_bam then num_rep_bam
        else length(nodup_bams)
    Int num_rep_ta = if length(tas)<num_rep_nodup_bam then num_rep_nodup_bam
        else length(tas)
    Int num_rep_peak = if length(peaks)<num_rep_ta then num_rep_ta
        else length(peaks)
    Int num_rep = num_rep_peak

    # temporary variables to get number of controls
    Int num_ctl_fastq = length(ctl_fastqs_R1)
    Int num_ctl_bam = if length(ctl_bams)<num_ctl_fastq then num_ctl_fastq
        else length(ctl_bams)
    Int num_ctl_nodup_bam = if length(ctl_nodup_bams)<num_ctl_bam then num_ctl_bam
        else length(ctl_nodup_bams)
    Int num_ctl_ta = if length(ctl_tas)<num_ctl_nodup_bam then num_ctl_nodup_bam
        else length(ctl_tas)
    Int num_ctl = num_ctl_ta

    # sanity check for inputs
    if ( num_rep == 0 && num_ctl == 0 ) {
        call raise_exception as error_input_data { input:
            msg = 'No FASTQ/BAM/TAG-ALIGN/PEAK defined in your input JSON. Check if your FASTQs are defined as "chip.fastqs_repX_RY". DO NOT MISS suffix _R1 even for single ended FASTQ.'
        }
    }
    if ( !align_only_ && peak_caller_ == 'spp' && num_ctl == 0 ) {
        call raise_exception as error_control_required { input:
            msg = 'SPP requires control inputs. Define control input files ("chip.ctl_*") in an input JSON file.'
        }
    }
    if ( (num_rep_fastq > 0 || num_ctl_fastq > 0) && aligner_ != 'bwa' && aligner_ != 'bowtie2' && aligner_ != 'custom' ) {
        call raise_exception as error_wrong_aligner { input:
            msg = 'Choose chip.aligner to align your fastqs. Choices: bwa, bowtie2, custom.'
        }
    }
    if ( aligner_ != 'bwa' && use_bwa_mem_for_pe ) {
        call raise_exception as error_use_bwa_mem_for_non_bwa { input:
            msg = 'To use chip.use_bwa_mem_for_pe, choose bwa for chip.aligner.'
        }
    }
    if ( aligner_ != 'bowtie2' && use_bowtie2_local_mode ) {
        call raise_exception as error_use_bowtie2_local_mode_for_non_bowtie2 { input:
            msg = 'To use chip.use_bowtie2_local_mode, choose bowtie2 for chip.aligner.'
        }
    }
    if ( aligner_ == 'custom' && ( !defined(custom_align_py) || !defined(custom_aligner_idx_tar) ) ) {
        call raise_exception as error_custom_aligner { input:
            msg = 'To use a custom aligner, define chip.custom_align_py and chip.custom_aligner_idx_tar.'
        }
    }

    if ( ( ctl_depth_limit > 0 || exp_ctl_depth_ratio_limit > 0 ) && num_ctl > 1 && length(ctl_paired_ends) > 1  ) {
        call raise_exception as error_subsample_pooled_control_with_mixed_endedness { input:
            msg = 'Cannot use automatic control subsampling ("chip.ctl_depth_limit">0 and "chip.exp_ctl_depth_limit">0) for ' +
                  'multiple controls with mixed endedness (e.g. SE ctl-rep1 and PE ctl-rep2). ' +
                  'Automatic control subsampling is enabled by default. ' +
                  'Disable automatic control subsampling by explicitly defining the above two parameters as 0 in your input JSON file. ' +
                  'You can still use manual control subsamping ("chip.ctl_subsample_reads">0) since it is done ' +
                  'for individual control\'s TAG-ALIGN output according to each control\'s endedness. '
        }
    }
    if ( pipeline_type == 'control' && num_ctl > 0 ) {
        call raise_exception as error_ctl_input_defined_in_control_mode { input:
            msg = 'In control mode (chip.pipeline_type: control), do not define ctl_* input variables. Define fastqs_repX_RY instead.'
        }
    }
    if ( pipeline_type == 'control' && num_rep_fastq == 0 ) {
        call raise_exception as error_ctl_fastq_input_required_for_control_mode { input:
            msg = 'Control mode (chip.pipeline_type: control) is for FASTQs only. Define FASTQs in fastqs_repX_RY. Pipeline will recognize them as control FASTQs.'
        }
    }

    # align each replicate
    scatter(i in range(num_rep)) {
        # to override endedness definition for individual replicate
        #     paired_end will override paired_ends[i]
        Boolean paired_end_ = if !defined(paired_end) && i<length(paired_ends) then paired_ends[i]
            else select_first([paired_end])

        Boolean has_input_of_align = i<length(fastqs_R1) && length(fastqs_R1[i])>0
        Boolean has_output_of_align = i<length(bams) && defined(bams[i])
        if ( has_input_of_align && !has_output_of_align ) {
            call align { input :
                fastqs_R1 = fastqs_R1[i],
                fastqs_R2 = fastqs_R2[i],
                crop_length = crop_length,
                crop_length_tol = crop_length_tol,
                trimmomatic_phred_score_format = trimmomatic_phred_score_format,

                aligner = aligner_,
                mito_chr_name = mito_chr_name_,
                custom_align_py = custom_align_py,
                idx_tar = if aligner=='bwa' then bwa_idx_tar_
                    else if aligner=='bowtie2' then bowtie2_idx_tar_
                    else custom_aligner_idx_tar,
                paired_end = paired_end_,
                use_bwa_mem_for_pe = use_bwa_mem_for_pe,
                bwa_mem_read_len_limit = bwa_mem_read_len_limit,
                use_bowtie2_local_mode = use_bowtie2_local_mode,
                ref_fa = ref_fa_,

                trimmomatic_java_heap = align_trimmomatic_java_heap,
                cpu = align_cpu,
                mem_factor = align_mem_factor_,
                time_hr = align_time_hr,
                disk_factor = align_disk_factor_,
            }
        }
        File? bam_ = if has_output_of_align then bams[i] else align.bam

        Boolean has_input_of_filter = has_output_of_align || defined(align.bam)
        Boolean has_output_of_filter = i<length(nodup_bams) && defined(nodup_bams[i])
        # skip if we already have output of this step
        if ( has_input_of_filter && !has_output_of_filter ) {
            call filter { input :
                bam = bam_,
                paired_end = paired_end_,
                ref_fa = ref_fa_,
                redact_nodup_bam = redact_nodup_bam,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = no_dup_removal,
                mito_chr_name = mito_chr_name_,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
            }
        }
        File? nodup_bam_ = if has_output_of_filter then nodup_bams[i] else filter.nodup_bam

        Boolean has_input_of_bam2ta = has_output_of_filter || defined(filter.nodup_bam)
        Boolean has_output_of_bam2ta = i<length(tas) && defined(tas[i])
        if ( has_input_of_bam2ta && !has_output_of_bam2ta ) {
            call bam2ta { input :
                bam = nodup_bam_,
                subsample = subsample_reads,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
            }
        }
        File? ta_ = if has_output_of_bam2ta then tas[i] else bam2ta.ta

        Boolean has_input_of_spr = has_output_of_bam2ta || defined(bam2ta.ta)
        if ( has_input_of_spr && !align_only_ && !true_rep_only ) {
            call spr { input :
                ta = ta_,
                paired_end = paired_end_,
                mem_factor = spr_mem_factor,
                disk_factor = spr_disk_factor,
            }
        }

        Boolean has_input_of_count_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
        if ( has_input_of_count_signal_track && enable_count_signal_track_ ) {
            # generate count signal track
            call count_signal_track { input :
                ta = ta_,
                chrsz = chrsz_,
            }
        }

        if ( enable_gc_bias_ && defined(nodup_bam_) && defined(ref_fa_) ) {
            call gc_bias { input :
                nodup_bam = nodup_bam_,
                ref_fa = ref_fa_,
                picard_java_heap = gc_bias_picard_java_heap,
            }
        }

        # special trimming/mapping for xcor (when starting from FASTQs)
        if ( has_input_of_align ) {
            call align as align_R1 { input :
                fastqs_R1 = fastqs_R1[i],
                fastqs_R2 = [],
                trim_bp = xcor_trim_bp,
                crop_length = 0,
                crop_length_tol = 0,
                trimmomatic_phred_score_format = trimmomatic_phred_score_format,

                aligner = aligner_,
                mito_chr_name = mito_chr_name_,
                custom_align_py = custom_align_py,
                idx_tar = if aligner=='bwa' then bwa_idx_tar_
                    else if aligner=='bowtie2' then bowtie2_idx_tar_
                    else custom_aligner_idx_tar,
                paired_end = false,
                use_bwa_mem_for_pe = use_bwa_mem_for_pe,
                bwa_mem_read_len_limit = bwa_mem_read_len_limit,
                use_bowtie2_local_mode = use_bowtie2_local_mode,
                ref_fa = ref_fa_,

                cpu = align_cpu,
                mem_factor = align_mem_factor_,
                time_hr = align_time_hr,
                disk_factor = align_disk_factor_,
            }
            # no bam deduping for xcor
            call filter as filter_R1 { input :
                bam = align_R1.bam,
                paired_end = false,
                redact_nodup_bam = false,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = true,
                mito_chr_name = mito_chr_name_,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
            }
            call bam2ta as bam2ta_no_dedup_R1 { input :
                bam = filter_R1.nodup_bam,  # it's named as nodup bam but it's not deduped but just filtered
                paired_end = false,
                subsample = 0,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
            }
        }

        # special trimming/mapping for xcor (when starting from BAMs)
        Boolean has_input_of_bam2ta_no_dedup = (has_output_of_align || defined(align.bam))
            && !defined(bam2ta_no_dedup_R1.ta)
        if ( has_input_of_bam2ta_no_dedup ) {
            call filter as filter_no_dedup { input :
                bam = bam_,
                paired_end = paired_end_,
                redact_nodup_bam = false,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = true,
                mito_chr_name = mito_chr_name_,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
            }
            call bam2ta as bam2ta_no_dedup { input :
                bam = filter_no_dedup.nodup_bam,  # output name is nodup but it's not deduped
                paired_end = paired_end_,
                subsample = 0,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
            }
        }

        # use trimmed/unfilitered R1 tagAlign for paired end dataset
        # if not starting from fastqs, keep using old method
        #  (mapping with both ends for tag-aligns to be used for xcor)
        # subsample tagalign (non-mito) and cross-correlation analysis
        File? ta_xcor = if defined(bam2ta_no_dedup_R1.ta) then bam2ta_no_dedup_R1.ta
            else if defined(bam2ta_no_dedup.ta) then bam2ta_no_dedup.ta
            else ta_
        Boolean paired_end_xcor = if defined(bam2ta_no_dedup_R1.ta) then false
            else paired_end_

        Boolean has_input_of_xcor = defined(ta_xcor)
        if ( has_input_of_xcor && enable_xcor_ ) {
            call xcor { input :
                ta = ta_xcor,
                paired_end = paired_end_xcor,
                subsample = xcor_subsample_reads,
                mito_chr_name = mito_chr_name_,
                chip_seq_type = pipeline_type,
                exclusion_range_min = xcor_exclusion_range_min,
                exclusion_range_max = xcor_exclusion_range_max,
                cpu = xcor_cpu,
                mem_factor = xcor_mem_factor,
                time_hr = xcor_time_hr,
                disk_factor = xcor_disk_factor,
            }
        }

        # before peak calling, get fragment length from xcor analysis or given input
        # if fraglen [] is defined in the input JSON, fraglen from xcor will be ignored
        Int? fraglen_ = if i<length(fraglen) && defined(fraglen[i]) then fraglen[i]
            else xcor.fraglen
    }

    # align each control
    scatter(i in range(num_ctl)) {
        # to override endedness definition for individual control
        #     ctl_paired_end will override ctl_paired_ends[i]
        Boolean ctl_paired_end_ = if !defined(ctl_paired_end) && i<length(ctl_paired_ends) then ctl_paired_ends[i]
            else select_first([ctl_paired_end, paired_end])

        Boolean has_input_of_align_ctl = i<length(ctl_fastqs_R1) && length(ctl_fastqs_R1[i])>0
        Boolean has_output_of_align_ctl = i<length(ctl_bams) && defined(ctl_bams[i])
        if ( has_input_of_align_ctl && !has_output_of_align_ctl ) {
            call align as align_ctl { input :
                fastqs_R1 = ctl_fastqs_R1[i],
                fastqs_R2 = ctl_fastqs_R2[i],
                crop_length = crop_length,
                crop_length_tol = crop_length_tol,
                trimmomatic_phred_score_format = trimmomatic_phred_score_format,

                aligner = aligner_,
                mito_chr_name = mito_chr_name_,
                custom_align_py = custom_align_py,
                idx_tar = if aligner=='bwa' then bwa_idx_tar_
                    else if aligner=='bowtie2' then bowtie2_idx_tar_
                    else custom_aligner_idx_tar,
                paired_end = ctl_paired_end_,
                use_bwa_mem_for_pe = use_bwa_mem_for_pe,
                bwa_mem_read_len_limit = bwa_mem_read_len_limit,
                use_bowtie2_local_mode = use_bowtie2_local_mode,
                ref_fa = ref_fa_,

                trimmomatic_java_heap = align_trimmomatic_java_heap,
                cpu = align_cpu,
                mem_factor = align_mem_factor_,
                time_hr = align_time_hr,
                disk_factor = align_disk_factor_,
            }
        }
        File? ctl_bam_ = if has_output_of_align_ctl then ctl_bams[i] else align_ctl.bam

        Boolean has_input_of_filter_ctl = has_output_of_align_ctl || defined(align_ctl.bam)
        Boolean has_output_of_filter_ctl = i<length(ctl_nodup_bams) && defined(ctl_nodup_bams[i])
        # skip if we already have output of this step
        if ( has_input_of_filter_ctl && !has_output_of_filter_ctl ) {
            call filter as filter_ctl { input :
                bam = ctl_bam_,
                paired_end = ctl_paired_end_,
                ref_fa = ref_fa_,
                redact_nodup_bam = redact_nodup_bam,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = no_dup_removal,
                mito_chr_name = mito_chr_name_,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
            }
        }
        File? ctl_nodup_bam_ = if has_output_of_filter_ctl then ctl_nodup_bams[i] else filter_ctl.nodup_bam

        Boolean has_input_of_bam2ta_ctl = has_output_of_filter_ctl || defined(filter_ctl.nodup_bam)
        Boolean has_output_of_bam2ta_ctl = i<length(ctl_tas) && defined(ctl_tas[i])
        if ( has_input_of_bam2ta_ctl && !has_output_of_bam2ta_ctl ) {
            call bam2ta as bam2ta_ctl { input :
                bam = ctl_nodup_bam_,
                subsample = subsample_reads,
                paired_end = ctl_paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
            }
        }
        File? ctl_ta_ = if has_output_of_bam2ta_ctl then ctl_tas[i] else bam2ta_ctl.ta
    }

    # if there are TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta = length(select_all(ta_))==num_rep
    if ( has_all_inputs_of_pool_ta && num_rep>1 ) {
        # pool tagaligns from true replicates
        call pool_ta { input :
            tas = ta_,
            prefix = 'rep',
        }
    }

    # if there are pr1 TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_pr1 = length(select_all(spr.ta_pr1))==num_rep
    if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only_ && !true_rep_only ) {
        # pool tagaligns from pseudo replicate 1
        call pool_ta as pool_ta_pr1 { input :
            tas = spr.ta_pr1,
            prefix = 'rep-pr1',
        }
    }

    # if there are pr2 TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_pr2 = length(select_all(spr.ta_pr2))==num_rep
    if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only_ && !true_rep_only ) {
        # pool tagaligns from pseudo replicate 2
        call pool_ta as pool_ta_pr2 { input :
            tas = spr.ta_pr2,
            prefix = 'rep-pr2',
        }
    }

    # if there are CTL TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_ctl = length(select_all(ctl_ta_))==num_ctl
    if ( has_all_inputs_of_pool_ta_ctl && num_ctl>1 ) {
        # pool tagaligns from true replicates
        call pool_ta as pool_ta_ctl { input :
            tas = ctl_ta_,
            prefix = 'ctl',
        }
    }

    Boolean has_input_of_count_signal_track_pooled = defined(pool_ta.ta_pooled)
    if ( has_input_of_count_signal_track_pooled && enable_count_signal_track_ && num_rep>1 ) {
        call count_signal_track as count_signal_track_pooled { input :
            ta = pool_ta.ta_pooled,
            chrsz = chrsz_,
        }
    }

    Boolean has_input_of_jsd = defined(blacklist_) && length(select_all(nodup_bam_))==num_rep
    if ( has_input_of_jsd && num_rep > 0 && enable_jsd_ ) {
        # fingerprint and JS-distance plot
        call jsd { input :
            nodup_bams = nodup_bam_,
            ctl_bams = ctl_nodup_bam_, # use first control only
            blacklist = blacklist_,
            mapq_thresh = mapq_thresh_,

            cpu = jsd_cpu,
            mem_factor = jsd_mem_factor,
            time_hr = jsd_time_hr,
            disk_factor = jsd_disk_factor,
        }
    }

    Boolean has_all_input_of_choose_ctl = length(select_all(ta_))==num_rep
        && length(select_all(ctl_ta_))==num_ctl && num_ctl > 0
    if ( has_all_input_of_choose_ctl && !align_only_ ) {
        # choose appropriate control for each exp IP replicate
        # outputs:
        #     choose_ctl.idx : control replicate index for each exp replicate 
        #                    -1 means pooled ctl replicate
        call choose_ctl { input:
            tas = ta_,
            ctl_tas = ctl_ta_,
            ta_pooled = pool_ta.ta_pooled,
            ctl_ta_pooled = pool_ta_ctl.ta_pooled,
            always_use_pooled_ctl = always_use_pooled_ctl,
            ctl_depth_ratio = ctl_depth_ratio,
            ctl_depth_limit = ctl_depth_limit,
            exp_ctl_depth_ratio_limit = exp_ctl_depth_ratio_limit,
        }
    }

    scatter(i in range(num_rep)) {
        # make control ta array [[1,2,3,4]] -> [[1],[2],[3],[4]]
        # chosen_ctl_ta_id
        #     >=0: control TA index (this means that control TA with this index exists)
        #     -1: use pooled control
        #    -2: there is no control
        Int chosen_ctl_ta_id = if has_all_input_of_choose_ctl && !align_only_ then
            select_first([choose_ctl.chosen_ctl_ta_ids])[i] else -2
        Int chosen_ctl_ta_subsample = if has_all_input_of_choose_ctl && !align_only_ then
            select_first([choose_ctl.chosen_ctl_ta_subsample])[i] else 0
        Boolean chosen_ctl_paired_end = if chosen_ctl_ta_id == -2 then false
            else if chosen_ctl_ta_id == -1 then ctl_paired_end_[0]
            else ctl_paired_end_[chosen_ctl_ta_id]

        if ( chosen_ctl_ta_id > -2 && chosen_ctl_ta_subsample > 0 ) {
            call subsample_ctl { input:
                ta = if chosen_ctl_ta_id == -1 then pool_ta_ctl.ta_pooled
                     else ctl_ta_[ chosen_ctl_ta_id ],
                subsample = chosen_ctl_ta_subsample,
                paired_end = chosen_ctl_paired_end,
                mem_factor = subsample_ctl_mem_factor,
                disk_factor = subsample_ctl_disk_factor,
            }
        }
        Array[File] chosen_ctl_tas = if chosen_ctl_ta_id <= -2 then []
            else if chosen_ctl_ta_subsample > 0 then [ select_first([subsample_ctl.ta_subsampled]) ]
            else if chosen_ctl_ta_id == -1 then [ select_first([pool_ta_ctl.ta_pooled]) ]
            else [ select_first([ctl_ta_[ chosen_ctl_ta_id ]]) ]
    }
    Int chosen_ctl_ta_pooled_subsample = if has_all_input_of_choose_ctl && !align_only_ then
        select_first([choose_ctl.chosen_ctl_ta_subsample_pooled]) else 0

    # workaround for dx error (Unsupported combination: womType: Int womValue: ([225], Array[Int]))
    Array[Int] fraglen_tmp = select_all(fraglen_)

    # we have all tas and ctl_tas (optional for histone chipseq) ready, let's call peaks
    scatter(i in range(num_rep)) {
        Boolean has_input_of_call_peak = defined(ta_[i])
        Boolean has_output_of_call_peak = i<length(peaks) && defined(peaks[i])
        if ( has_input_of_call_peak && !has_output_of_call_peak && !align_only_ ) {
            call call_peak { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                tas = flatten([[ta_[i]], chosen_ctl_tas[i]]),
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                fdr_thresh = fdr_thresh,
                fraglen = fraglen_tmp[i],
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor_,
                disk_factor = call_peak_disk_factor_,
                time_hr = call_peak_time_hr,
            }
        }
        File? peak_ = if has_output_of_call_peak then peaks[i]
            else call_peak.peak

        # signal track
        if ( has_input_of_call_peak && !align_only_ ) {
            call macs2_signal_track { input :
                tas = flatten([[ta_[i]], chosen_ctl_tas[i]]),
                gensz = gensz_,
                chrsz = chrsz_,
                pval_thresh = pval_thresh,
                fraglen = fraglen_tmp[i],

                mem_factor = macs2_signal_track_mem_factor,
                disk_factor = macs2_signal_track_disk_factor,
                time_hr = macs2_signal_track_time_hr,
            }
        }

        # call peaks on 1st pseudo replicated tagalign
        Boolean has_input_of_call_peak_pr1 = defined(spr.ta_pr1[i])
        Boolean has_output_of_call_peak_pr1 = i<length(peaks_pr1) && defined(peaks_pr1[i])
        if ( has_input_of_call_peak_pr1 && !has_output_of_call_peak_pr1 && !true_rep_only ) {
            call call_peak as call_peak_pr1 { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                tas = flatten([[spr.ta_pr1[i]], chosen_ctl_tas[i]]),
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                fdr_thresh = fdr_thresh,
                fraglen = fraglen_tmp[i],
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
    
                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor_,
                disk_factor = call_peak_disk_factor_,
                time_hr = call_peak_time_hr,
            }
        }
        File? peak_pr1_ = if has_output_of_call_peak_pr1 then peaks_pr1[i]
            else call_peak_pr1.peak

        # call peaks on 2nd pseudo replicated tagalign
        Boolean has_input_of_call_peak_pr2 = defined(spr.ta_pr2[i])
        Boolean has_output_of_call_peak_pr2 = i<length(peaks_pr2) && defined(peaks_pr2[i])
        if ( has_input_of_call_peak_pr2 && !has_output_of_call_peak_pr2 && !true_rep_only ) {
            call call_peak as call_peak_pr2 { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                tas = flatten([[spr.ta_pr2[i]], chosen_ctl_tas[i]]),
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                fdr_thresh = fdr_thresh,
                fraglen = fraglen_tmp[i],
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor_,
                disk_factor = call_peak_disk_factor_,
                time_hr = call_peak_time_hr,
            }
        }
        File? peak_pr2_ = if has_output_of_call_peak_pr2 then peaks_pr2[i]
            else call_peak_pr2.peak
    }

    # if ( !align_only_ && num_rep > 1 ) {
    # rounded mean of fragment length, which will be used for 
    #  1) calling peaks for pooled true/pseudo replicates
    #  2) calculating FRiP
    call rounded_mean as fraglen_mean { input :
        ints = fraglen_tmp,
    }
    # }

    if ( has_all_input_of_choose_ctl && !align_only_ && chosen_ctl_ta_pooled_subsample > 0 ) {
        call subsample_ctl as subsample_ctl_pooled { input:
            ta = if num_ctl < 2 then ctl_ta_[0]
                 else pool_ta_ctl.ta_pooled,
            subsample = chosen_ctl_ta_pooled_subsample,
            paired_end = ctl_paired_end_[0],
            mem_factor = subsample_ctl_mem_factor,
            disk_factor = subsample_ctl_disk_factor,
        }
    }
    # actually not an array
    Array[File?] chosen_ctl_ta_pooled = if !has_all_input_of_choose_ctl || align_only_ then []
        else if chosen_ctl_ta_pooled_subsample > 0 then [ subsample_ctl_pooled.ta_subsampled ]
        else if num_ctl < 2 then [ ctl_ta_[0] ]
        else [ pool_ta_ctl.ta_pooled ]

    Boolean has_input_of_call_peak_pooled = defined(pool_ta.ta_pooled)
    Boolean has_output_of_call_peak_pooled = defined(peak_pooled)
    if ( has_input_of_call_peak_pooled && !has_output_of_call_peak_pooled && !align_only_ && num_rep>1 ) {
        # call peaks on pooled replicate
        # always call peaks for pooled replicate to get signal tracks
        call call_peak as call_peak_pooled { input :
            peak_caller = peak_caller_,
            peak_type = peak_type_,
            tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
            gensz = gensz_,
            chrsz = chrsz_,
            cap_num_peak = cap_num_peak_,
            pval_thresh = pval_thresh,
            fdr_thresh = fdr_thresh,
            fraglen = fraglen_mean.rounded_mean,
            blacklist = blacklist_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

            cpu = call_peak_cpu,
            mem_factor = call_peak_mem_factor_,
            disk_factor = call_peak_disk_factor_,
            time_hr = call_peak_time_hr,
        }
    }
    File? peak_pooled_ = if has_output_of_call_peak_pooled then peak_pooled
        else call_peak_pooled.peak    

    # macs2 signal track for pooled rep
    if ( has_input_of_call_peak_pooled && !align_only_ && num_rep>1 ) {
        call macs2_signal_track as macs2_signal_track_pooled { input :
            tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
            gensz = gensz_,
            chrsz = chrsz_,
            pval_thresh = pval_thresh,
            fraglen = fraglen_mean.rounded_mean,

            mem_factor = macs2_signal_track_mem_factor,
            disk_factor = macs2_signal_track_disk_factor,
            time_hr = macs2_signal_track_time_hr,
        }
    }

    Boolean has_input_of_call_peak_ppr1 = defined(pool_ta_pr1.ta_pooled)
    Boolean has_output_of_call_peak_ppr1 = defined(peak_ppr1)
    if ( has_input_of_call_peak_ppr1 && !has_output_of_call_peak_ppr1 && !align_only_ && !true_rep_only && num_rep>1 ) {
        # call peaks on 1st pooled pseudo replicates
        call call_peak as call_peak_ppr1 { input :
            peak_caller = peak_caller_,
            peak_type = peak_type_,
            tas = flatten([select_all([pool_ta_pr1.ta_pooled]), chosen_ctl_ta_pooled]),
            gensz = gensz_,
            chrsz = chrsz_,
            cap_num_peak = cap_num_peak_,
            pval_thresh = pval_thresh,
            fdr_thresh = fdr_thresh,
            fraglen = fraglen_mean.rounded_mean,
            blacklist = blacklist_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

            cpu = call_peak_cpu,
            mem_factor = call_peak_mem_factor_,
            disk_factor = call_peak_disk_factor_,
            time_hr = call_peak_time_hr,
        }
    }
    File? peak_ppr1_ = if has_output_of_call_peak_ppr1 then peak_ppr1
        else call_peak_ppr1.peak

    Boolean has_input_of_call_peak_ppr2 = defined(pool_ta_pr2.ta_pooled)
    Boolean has_output_of_call_peak_ppr2 = defined(peak_ppr2)
    if ( has_input_of_call_peak_ppr2 && !has_output_of_call_peak_ppr2 && !align_only_ && !true_rep_only && num_rep>1 ) {
        # call peaks on 2nd pooled pseudo replicates
        call call_peak as call_peak_ppr2 { input :
            peak_caller = peak_caller_,
            peak_type = peak_type_,
            tas = flatten([select_all([pool_ta_pr2.ta_pooled]), chosen_ctl_ta_pooled]),
            gensz = gensz_,
            chrsz = chrsz_,
            cap_num_peak = cap_num_peak_,
            pval_thresh = pval_thresh,
            fdr_thresh = fdr_thresh,
            fraglen = fraglen_mean.rounded_mean,
            blacklist = blacklist_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

            cpu = call_peak_cpu,
            mem_factor = call_peak_mem_factor_,
            disk_factor = call_peak_disk_factor_,
            time_hr = call_peak_time_hr,
        }
    }
    File? peak_ppr2_ = if has_output_of_call_peak_ppr2 then peak_ppr2
        else call_peak_ppr2.peak

    # do IDR/overlap on all pairs of two replicates (i,j)
    #     where i and j are zero-based indices and 0 <= i < j < num_rep
    scatter( pair in cross(range(num_rep),range(num_rep)) ) {
        # pair.left = 0-based index of 1st replicate
        # pair.right = 0-based index of 2nd replicate
        File? peak1_ = peak_[pair.left]
        File? peak2_ = peak_[pair.right]
        if ( !align_only_ && pair.left<pair.right ) {
            # Naive overlap on every pair of true replicates
            call overlap { input :
                prefix = 'rep'+(pair.left+1)+'_vs_rep'+(pair.right+1),
                peak1 = peak1_,
                peak2 = peak2_,
                peak_pooled = peak_pooled_,
                fraglen = fraglen_mean.rounded_mean,
                peak_type = peak_type_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = pool_ta.ta_pooled,
            }
        }
        if ( enable_idr && !align_only_ && pair.left<pair.right ) {
            # IDR on every pair of true replicates
            call idr { input :
                prefix = 'rep'+(pair.left+1)+'_vs_rep'+(pair.right+1),
                peak1 = peak1_,
                peak2 = peak2_,
                peak_pooled = peak_pooled_,
                fraglen = fraglen_mean.rounded_mean,
                idr_thresh = idr_thresh,
                peak_type = peak_type_,
                rank = idr_rank_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = pool_ta.ta_pooled,
            }
        }
    }

    # overlap on pseudo-replicates (pr1, pr2) for each true replicate
    if ( !align_only_ && !true_rep_only ) {
        scatter( i in range(num_rep) ) {
            call overlap as overlap_pr { input :
                prefix = 'rep'+(i+1)+'-pr1_vs_rep'+(i+1)+'-pr2',
                peak1 = peak_pr1_[i],
                peak2 = peak_pr2_[i],
                peak_pooled = peak_[i],
                fraglen = fraglen_[i],
                peak_type = peak_type_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = ta_[i],
            }
        }
    }

    if ( !align_only_ && !true_rep_only && enable_idr ) {
        scatter( i in range(num_rep) ) {
            # IDR on pseduo replicates
            call idr as idr_pr { input :
                prefix = 'rep'+(i+1)+'-pr1_vs_rep'+(i+1)+'-pr2',
                peak1 = peak_pr1_[i],
                peak2 = peak_pr2_[i],
                peak_pooled = peak_[i],
                fraglen = fraglen_[i],
                idr_thresh = idr_thresh,
                peak_type = peak_type_,
                rank = idr_rank_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = ta_[i],
            }
        }
    }

    if ( !align_only_ && !true_rep_only && num_rep > 1 ) {
        # Naive overlap on pooled pseudo replicates
        call overlap as overlap_ppr { input :
            prefix = 'pooled-pr1_vs_pooled-pr2',
            peak1 = peak_ppr1_,
            peak2 = peak_ppr2_,
            peak_pooled = peak_pooled_,
            peak_type = peak_type_,
            fraglen = fraglen_mean.rounded_mean,
            blacklist = blacklist_,
            chrsz = chrsz_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
            ta = pool_ta.ta_pooled,
        }
    }

    if ( !align_only_ && !true_rep_only && num_rep > 1 && enable_idr ) {
        # IDR on pooled pseduo replicates
        call idr as idr_ppr { input :
            prefix = 'pooled-pr1_vs_pooled-pr2',
            peak1 = peak_ppr1_,
            peak2 = peak_ppr2_,
            peak_pooled = peak_pooled_,
            idr_thresh = idr_thresh,
            peak_type = peak_type_,
            fraglen = fraglen_mean.rounded_mean,
            rank = idr_rank_,
            blacklist = blacklist_,
            chrsz = chrsz_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
            ta = pool_ta.ta_pooled,
        }
    }

    # reproducibility QC for overlap/IDR peaks
    if ( !align_only_ && !true_rep_only && num_rep > 0 ) {
        # reproducibility QC for overlapping peaks
        call reproducibility as reproducibility_overlap { input :
            prefix = 'overlap',
            peaks = select_all(overlap.bfilt_overlap_peak),
            peaks_pr = overlap_pr.bfilt_overlap_peak,
            peak_ppr = overlap_ppr.bfilt_overlap_peak,
            peak_type = peak_type_,
            chrsz = chrsz_,
        }
    }

    if ( !align_only_ && !true_rep_only && num_rep > 0 && enable_idr ) {
        # reproducibility QC for IDR peaks
        call reproducibility as reproducibility_idr { input :
            prefix = 'idr',
            peaks = select_all(idr.bfilt_idr_peak),
            peaks_pr = idr_pr.bfilt_idr_peak,
            peak_ppr = idr_ppr.bfilt_idr_peak,
            peak_type = peak_type_,
            chrsz = chrsz_,
        }
    }

    # Generate final QC report and JSON
    call qc_report { input :
        pipeline_ver = pipeline_ver,
        title = title,
        description = description,
        genome = genome_name_,
        paired_ends = paired_end_,
        ctl_paired_ends = ctl_paired_end_,
        pipeline_type = pipeline_type,
        aligner = aligner_,
        peak_caller = peak_caller_,
        cap_num_peak = cap_num_peak_,
        idr_thresh = idr_thresh,
        pval_thresh = pval_thresh,
        xcor_trim_bp = xcor_trim_bp,
        xcor_subsample_reads = xcor_subsample_reads,

        samstat_qcs = select_all(align.samstat_qc),
        nodup_samstat_qcs = select_all(filter.samstat_qc),
        dup_qcs = select_all(filter.dup_qc),
        lib_complexity_qcs = select_all(filter.lib_complexity_qc),
        xcor_plots = select_all(xcor.plot_png),
        xcor_scores = select_all(xcor.score),

        ctl_samstat_qcs = select_all(align_ctl.samstat_qc),
        ctl_nodup_samstat_qcs = select_all(filter_ctl.samstat_qc),
        ctl_dup_qcs = select_all(filter_ctl.dup_qc),
        ctl_lib_complexity_qcs = select_all(filter_ctl.lib_complexity_qc),

        jsd_plot = jsd.plot,
        jsd_qcs = jsd.jsd_qcs,

        frip_qcs = select_all(call_peak.frip_qc),
        frip_qcs_pr1 = select_all(call_peak_pr1.frip_qc),
        frip_qcs_pr2 = select_all(call_peak_pr2.frip_qc),
        frip_qc_pooled = call_peak_pooled.frip_qc,
        frip_qc_ppr1 = call_peak_ppr1.frip_qc,
        frip_qc_ppr2 = call_peak_ppr2.frip_qc,

        idr_plots = select_all(idr.idr_plot),
        idr_plots_pr = idr_pr.idr_plot,
        idr_plot_ppr = idr_ppr.idr_plot,
        frip_idr_qcs = select_all(idr.frip_qc),
        frip_idr_qcs_pr = idr_pr.frip_qc,
        frip_idr_qc_ppr = idr_ppr.frip_qc,
        frip_overlap_qcs = select_all(overlap.frip_qc),
        frip_overlap_qcs_pr = overlap_pr.frip_qc,
        frip_overlap_qc_ppr = overlap_ppr.frip_qc,
        idr_reproducibility_qc = reproducibility_idr.reproducibility_qc,
        overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,

        gc_plots = select_all(gc_bias.gc_plot),

        peak_region_size_qcs = select_all(call_peak.peak_region_size_qc),
        peak_region_size_plots = select_all(call_peak.peak_region_size_plot),
        num_peak_qcs = select_all(call_peak.num_peak_qc),

        idr_opt_peak_region_size_qc = reproducibility_idr.peak_region_size_qc,
        idr_opt_peak_region_size_plot = reproducibility_overlap.peak_region_size_plot,
        idr_opt_num_peak_qc = reproducibility_idr.num_peak_qc,

        overlap_opt_peak_region_size_qc = reproducibility_overlap.peak_region_size_qc,
        overlap_opt_peak_region_size_plot = reproducibility_overlap.peak_region_size_plot,
        overlap_opt_num_peak_qc = reproducibility_overlap.num_peak_qc,
    }

    output {
        File report = qc_report.report
        File qc_json = qc_report.qc_json
        Boolean qc_json_ref_match = qc_report.qc_json_ref_match
    }
}

task align {
    input {
        Array[File] fastqs_R1         # [merge_id]
        Array[File] fastqs_R2
        File? ref_fa
        Int? trim_bp            # this is for R1 only
        Int crop_length
        Int crop_length_tol
        String? trimmomatic_phred_score_format

        String aligner

        String mito_chr_name
        Int? multimapping
        File? custom_align_py
        File? idx_tar            # reference index tar
        Boolean paired_end
        Boolean use_bwa_mem_for_pe
        Int bwa_mem_read_len_limit
        Boolean use_bowtie2_local_mode

        String? trimmomatic_java_heap
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(fastqs_R1, "G") + size(fastqs_R2, "G")
    Float mem_gb = 5.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    Float trimmomatic_java_heap_factor = 0.9
    Array[Array[File]] tmp_fastqs = if paired_end then transpose([fastqs_R1, fastqs_R2])
                else transpose([fastqs_R1])
    command {
        set -e

        # check if pipeline dependencies can be found
        if [[ -z "$(which encode_task_merge_fastq.py 2> /dev/null || true)" ]]
        then
          echo -e "\n* Error: pipeline dependencies not found." 1>&2
          echo 'Conda users: Did you activate Conda environment (conda activate encode-chip-seq-pipeline)?' 1>&2
          echo '    Or did you install Conda and environment correctly (bash scripts/install_conda_env.sh)?' 1>&2
          echo 'GCP/AWS/Docker users: Did you add --docker flag to Caper command line arg?' 1>&2
          echo 'Singularity users: Did you add --singularity flag to Caper command line arg?' 1>&2
          echo -e "\n" 1>&2
          exit 3
        fi
        python3 $(which encode_task_merge_fastq.py) \
            ${write_tsv(tmp_fastqs)} \
            ${if paired_end then '--paired-end' else ''} \
            ${'--nth ' + cpu}

        if [ -z '${trim_bp}' ]; then
            SUFFIX=
        else
            SUFFIX=_trimmed
            python3 $(which encode_task_trim_fastq.py) \
                R1/*.fastq.gz \
                --trim-bp ${trim_bp} \
                --out-dir R1$SUFFIX
            if [ '${paired_end}' == 'true' ]; then
                python3 $(which encode_task_trim_fastq.py) \
                    R2/*.fastq.gz \
                    --trim-bp ${trim_bp} \
                    --out-dir R2$SUFFIX
            fi
        fi
        if [ '${crop_length}' == '0' ]; then
            SUFFIX=$SUFFIX
        else
            NEW_SUFFIX="$SUFFIX"_cropped
            python3 $(which encode_task_trimmomatic.py) \
                --fastq1 R1$SUFFIX/*.fastq.gz \
                ${if paired_end then '--fastq2 R2$SUFFIX/*.fastq.gz' else ''} \
                ${if paired_end then '--paired-end' else ''} \
                --crop-length ${crop_length} \
                --crop-length-tol "${crop_length_tol}" \
                ${'--phred-score-format ' + trimmomatic_phred_score_format } \
                --out-dir-R1 R1$NEW_SUFFIX \
                ${if paired_end then '--out-dir-R2 R2$NEW_SUFFIX' else ''} \
                ${'--trimmomatic-java-heap ' + if defined(trimmomatic_java_heap) then trimmomatic_java_heap else (round(mem_gb * trimmomatic_java_heap_factor) + 'G')} \
                ${'--nth ' + cpu}
            SUFFIX=$NEW_SUFFIX
        fi

        if [ '${aligner}' == 'bwa' ]; then
            python3 $(which encode_task_bwa.py) \
                ${idx_tar} \
                R1$SUFFIX/*.fastq.gz \
                ${if paired_end then 'R2$SUFFIX/*.fastq.gz' else ''} \
                ${if paired_end then '--paired-end' else ''} \
                ${if use_bwa_mem_for_pe then '--use-bwa-mem-for-pe' else ''} \
                ${'--bwa-mem-read-len-limit ' + bwa_mem_read_len_limit} \
                ${'--mem-gb ' + samtools_mem_gb} \
                ${'--nth ' + cpu}

        elif [ '${aligner}' == 'bowtie2' ]; then
            python3 $(which encode_task_bowtie2.py) \
                ${idx_tar} \
                R1$SUFFIX/*.fastq.gz \
                ${if paired_end then 'R2$SUFFIX/*.fastq.gz' else ''} \
                ${'--multimapping ' + multimapping} \
                ${if paired_end then '--paired-end' else ''} \
                ${if use_bowtie2_local_mode then '--local' else ''} \
                ${'--mem-gb ' + samtools_mem_gb} \
                ${'--nth ' + cpu}
        else
            python3 ${custom_align_py} \
                ${idx_tar} \
                R1$SUFFIX/*.fastq.gz \
                ${if paired_end then 'R2$SUFFIX/*.fastq.gz' else ''} \
                ${if paired_end then '--paired-end' else ''} \
                ${'--mem-gb ' + samtools_mem_gb} \
                ${'--nth ' + cpu}
        fi 

        python3 $(which encode_task_post_align.py) \
            R1$SUFFIX/*.fastq.gz $(ls *.bam) \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu}
        rm -rf R1 R2 R1$SUFFIX R2$SUFFIX
    }
    output {
        File bam = glob('*.bam')[0]
        File bai = glob('*.bai')[0]
        File samstat_qc = glob('*.samstats.qc')[0]
        File read_len_log = glob('*.read_length.txt')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0
    }
}

task filter {
    input {
        File? bam
        Boolean paired_end
        File? ref_fa
        Boolean redact_nodup_bam
        String dup_marker             # picard.jar MarkDuplicates (picard) or 
                                    # sambamba markdup (sambamba)
        Int mapq_thresh                # threshold for low MAPQ reads removal
        Array[String] filter_chrs     # chrs to be removed from final (nodup/filt) BAM
        File chrsz                    # 2-col chromosome sizes file
        Boolean no_dup_removal         # no dupe reads removal when filtering BAM
        String mito_chr_name

        Int cpu
        Float mem_factor
        String? picard_java_heap
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(bam, "G")
    Float picard_java_heap_factor = 0.9
    Float mem_gb = 6.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_filter.py) \
            ${bam} \
            ${if paired_end then '--paired-end' else ''} \
            --multimapping 0 \
            ${'--dup-marker ' + dup_marker} \
            ${'--mapq-thresh ' + mapq_thresh} \
            --filter-chrs ${sep=' ' filter_chrs} \
            ${'--chrsz ' + chrsz} \
            ${if no_dup_removal then '--no-dup-removal' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}

        if [ '${redact_nodup_bam}' == 'true' ]; then
            python3 $(which encode_task_bam_to_pbam.py) \
                $(ls *.bam) \
                ${'--ref-fa ' + ref_fa} \
                '--delete-original-bam'
        fi
    }
    output {
        File nodup_bam = glob('*.bam')[0]
        File nodup_bai = glob('*.bai')[0]
        File samstat_qc = glob('*.samstats.qc')[0]
        File dup_qc = glob('*.dup.qc')[0]
        File lib_complexity_qc = glob('*.lib_complexity.qc')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
    }
}

task bam2ta {
    input {
        File? bam
        Boolean paired_end
        String mito_chr_name         # mito chromosome name
        Int subsample                 # number of reads to subsample TAGALIGN
                                    # this affects all downstream analysis
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_bam2ta.py) \
            ${bam} \
            --disable-tn5-shift \
            ${if paired_end then '--paired-end' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--subsample ' + subsample} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu}
    }
    output {
        File ta = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
    }
}

task spr {
    input {
        File? ta
        Boolean paired_end

        Float mem_factor
        Float disk_factor
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_spr.py) \
            ${ta} \
            ${if paired_end then '--paired-end' else ''}
    }
    output {
        File ta_pr1 = glob('*.pr1.tagAlign.gz')[0]
        File ta_pr2 = glob('*.pr2.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 1
        disks : 'local-disk ${disk_gb} SSD'
    }
}

task pool_ta {
    input {
        Array[File?] tas
        Int? col             # number of columns in pooled TA
        String? prefix         # basename prefix
    }

    command {
        set -e
        python3 $(which encode_task_pool_ta.py) \
            ${sep=' ' select_all(tas)} \
            ${'--prefix ' + prefix} \
            ${'--col ' + col}
    }
    output {
        File ta_pooled = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 1
        disks : 'local-disk 100 SSD'
    }
}

task xcor {
    input {
        File? ta
        Boolean paired_end
        String mito_chr_name
        Int subsample  # number of reads to subsample TAGALIGN
                        # this will be used for xcor only
                        # will not affect any downstream analysis
        String? chip_seq_type
        Int? exclusion_range_min
        Int? exclusion_range_max

        Int cpu
        Float mem_factor    
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 8.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_xcor.py) \
            ${ta} \
            ${if paired_end then '--paired-end' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--subsample ' + subsample} \
            ${'--chip-seq-type ' + chip_seq_type} \
            ${'--exclusion-range-min ' + exclusion_range_min} \
            ${'--exclusion-range-max ' + exclusion_range_max} \
            ${'--subsample ' + subsample} \
            ${'--nth ' + cpu}
    }
    output {
        File plot_pdf = glob('*.cc.plot.pdf')[0]
        File plot_png = glob('*.cc.plot.png')[0]
        File score = glob('*.cc.qc')[0]
        File fraglen_log = glob('*.cc.fraglen.txt')[0]
        Int fraglen = read_int(fraglen_log)
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
    }
}

task jsd {
    input {
        Array[File?] nodup_bams
        Array[File?] ctl_bams
        File? blacklist
        Int mapq_thresh

        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(nodup_bams, "G") + size(ctl_bams, "G")
    Float mem_gb = 5.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_jsd.py) \
            ${sep=' ' select_all(nodup_bams)} \
            ${if length(ctl_bams)>0 then '--ctl-bam '+ select_first(ctl_bams) else ''} \
            ${'--mapq-thresh '+ mapq_thresh} \
            ${'--blacklist '+ blacklist} \
            ${'--nth ' + cpu}
    }
    output {
        File plot = glob('*.png')[0]
        Array[File] jsd_qcs = glob('*.jsd.qc')
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
    }
}

task choose_ctl {
    input {
        Array[File?] tas
        Array[File?] ctl_tas
        File? ta_pooled
        File? ctl_ta_pooled
        Boolean always_use_pooled_ctl # always use pooled control for all exp rep.
        Float ctl_depth_ratio         # if ratio between controls is higher than this
                                    # then always use pooled control for all exp rep.
        Int ctl_depth_limit
        Float exp_ctl_depth_ratio_limit
    }

    command {
        set -e
        python3 $(which encode_task_choose_ctl.py) \
            --tas ${sep=' ' select_all(tas)} \
            --ctl-tas ${sep=' ' select_all(ctl_tas)} \
            ${'--ta-pooled ' + ta_pooled} \
            ${'--ctl-ta-pooled ' + ctl_ta_pooled} \
            ${if always_use_pooled_ctl then '--always-use-pooled-ctl' else ''} \
            ${'--ctl-depth-ratio ' + ctl_depth_ratio} \
            ${'--ctl-depth-limit ' + ctl_depth_limit} \
            ${'--exp-ctl-depth-ratio-limit ' + exp_ctl_depth_ratio_limit}
    }
    output {
        File chosen_ctl_id_tsv = glob('chosen_ctl.tsv')[0]
        File chosen_ctl_subsample_tsv = glob('chosen_ctl_subsample.tsv')[0]
        File chosen_ctl_subsample_pooled_txt = glob('chosen_ctl_subsample_pooled.txt')[0]
        Array[Int] chosen_ctl_ta_ids = read_lines(chosen_ctl_id_tsv)
        Array[Int] chosen_ctl_ta_subsample = read_lines(chosen_ctl_subsample_tsv)
        Int chosen_ctl_ta_subsample_pooled = read_int(chosen_ctl_subsample_pooled_txt)
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 1
        disks : 'local-disk 50 SSD'
    }
}

task count_signal_track {
    input {
        File? ta             # tag-align
        File chrsz            # 2-col chromosome sizes file
    }

    command {
        set -e
        python3 $(which encode_task_count_signal_track.py) \
            ${ta} \
            ${'--chrsz ' + chrsz}
    }
    output {
        File pos_bw = glob('*.positive.bigwig')[0]
        File neg_bw = glob('*.negative.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 50 SSD'
    }
}

task subsample_ctl {
    input {
        File? ta
        Boolean paired_end
        Int subsample

        Float mem_factor
        Float disk_factor
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        python3 $(which encode_task_subsample_ctl.py) \
            ${ta} \
            ${'--subsample ' + subsample} \
            ${if paired_end then '--paired-end' else ''} \
    }
    output {
        File ta_subsampled = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 4
        disks : 'local-disk ${disk_gb} SSD'
    }
}

task call_peak {
    input {
        String peak_caller
        String peak_type
        Array[File?] tas    # [ta, control_ta]. control_ta is optional
        Int fraglen         # fragment length from xcor
        String gensz        # Genome size (sum of entries in 2nd column of 
                            # chr. sizes file, or hs for human, ms for mouse)
        File chrsz            # 2-col chromosome sizes file
        Int cap_num_peak    # cap number of raw peaks called from MACS2
        Float pval_thresh     # p.value threshold for MACS2
        Float? fdr_thresh     # FDR threshold for SPP

        File? blacklist     # blacklist BED to filter raw peaks
        String? regex_bfilt_peak_chr_name

        Int cpu    
        Float mem_factor
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(tas, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e

        if [ '${peak_caller}' == 'macs2' ]; then
            python3 $(which encode_task_macs2_chip.py) \
                ${sep=' ' select_all(tas)} \
                ${'--gensz '+ gensz} \
                ${'--chrsz ' + chrsz} \
                ${'--fraglen ' + fraglen} \
                ${'--cap-num-peak ' + cap_num_peak} \
                ${'--pval-thresh '+ pval_thresh}

        elif [ '${peak_caller}' == 'spp' ]; then
            python3 $(which encode_task_spp.py) \
                ${sep=' ' select_all(tas)} \
                ${'--chrsz ' + chrsz} \
                ${'--fraglen ' + fraglen} \
                ${'--cap-num-peak ' + cap_num_peak} \
                ${'--fdr-thresh '+ fdr_thresh} \
                ${'--nth ' + cpu}
        fi

        python3 $(which encode_task_post_call_peak_chip.py) \
            $(ls *Peak.gz) \
            ${'--ta ' + tas[0]} \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--chrsz ' + chrsz} \
            ${'--fraglen ' + fraglen} \
            ${'--peak-type ' + peak_type} \
            ${'--blacklist ' + blacklist}        
    }
    output {
        File peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        # generated by post_call_peak py
        File bfilt_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = glob('*.frip.qc')[0]
        File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
        File peak_region_size_plot = glob('*.peak_region_size.png')[0]
        File num_peak_qc = glob('*.num_peak.qc')[0]
    }
    runtime {
        cpu : if peak_caller == 'macs2' then 2 else cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0        
    }
}

task macs2_signal_track {
    input {
        Array[File?] tas    # [ta, control_ta]. control_ta is optional
        Int fraglen         # fragment length from xcor
        String gensz        # Genome size (sum of entries in 2nd column of 
                            # chr. sizes file, or hs for human, ms for mouse)
        File chrsz            # 2-col chromosome sizes file
        Float pval_thresh     # p.value threshold

        Float mem_factor
        Int time_hr
        Float disk_factor
    }
    Float input_file_size_gb = size(tas, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_macs2_signal_track_chip.py) \
            ${sep=' ' select_all(tas)} \
            ${'--gensz '+ gensz} \
            ${'--chrsz ' + chrsz} \
            ${'--fraglen ' + fraglen} \
            ${'--pval-thresh '+ pval_thresh}
    }
    output {
        File pval_bw = glob('*.pval.signal.bigwig')[0]
        File fc_bw = glob('*.fc.signal.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0
    }
}

task idr {
    input {
        String prefix         # prefix for IDR output file
        File? peak1
        File? peak2
        File? peak_pooled
        Float idr_thresh
        File? blacklist     # blacklist BED to filter raw peaks
        String regex_bfilt_peak_chr_name
        # parameters to compute FRiP
        File? ta            # to calculate FRiP
        Int? fraglen         # fragment length from xcor
        File chrsz            # 2-col chromosome sizes file
        String peak_type
        String rank
    }

    command {
        set -e
        ${if defined(ta) then '' else 'touch null.frip.qc'}
        touch null 
        python3 $(which encode_task_idr.py) \
            ${peak1} ${peak2} ${peak_pooled} \
            ${'--prefix ' + prefix} \
            ${'--idr-thresh ' + idr_thresh} \
            ${'--peak-type ' + peak_type} \
            --idr-rank ${rank} \
            ${'--fraglen ' + fraglen} \
            ${'--chrsz ' + chrsz} \
            ${'--blacklist '+ blacklist} \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--ta ' + ta}
    }
    output {
        File idr_peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        File bfilt_idr_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_idr_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_idr_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_idr_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_idr_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File idr_plot = glob('*.txt.png')[0]
        File idr_unthresholded_peak = glob('*.txt.gz')[0]
        File idr_log = glob('*.idr*.log')[0]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 1
        disks : 'local-disk 50 SSD'
    }    
}

task overlap {
    input {
        String prefix     # prefix for IDR output file
        File? peak1
        File? peak2
        File? peak_pooled
        File? blacklist # blacklist BED to filter raw peaks
        String regex_bfilt_peak_chr_name
        # parameters to compute FRiP
        File? ta        # to calculate FRiP
        Int? fraglen     # fragment length from xcor (for FRIP)
        File chrsz        # 2-col chromosome sizes file
        String peak_type
    }

    command {
        set -e
        ${if defined(ta) then '' else 'touch null.frip.qc'}
        touch null 
        python3 $(which encode_task_overlap.py) \
            ${peak1} ${peak2} ${peak_pooled} \
            ${'--prefix ' + prefix} \
            ${'--peak-type ' + peak_type} \
            ${'--fraglen ' + fraglen} \
            ${'--chrsz ' + chrsz} \
            ${'--blacklist '+ blacklist} \
            --nonamecheck \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--ta ' + ta}
    }
    output {
        File overlap_peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        File bfilt_overlap_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_overlap_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_overlap_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_overlap_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_overlap_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 1
        disks : 'local-disk 50 SSD'
    }    
}

task reproducibility {
    input {
        String prefix
        Array[File] peaks # peak files from pair of true replicates
                            # in a sorted order. for example of 4 replicates,
                            # 1,2 1,3 1,4 2,3 2,4 3,4.
                            # x,y means peak file from rep-x vs rep-y
        Array[File]? peaks_pr    # peak files from pseudo replicates
        File? peak_ppr            # Peak file from pooled pseudo replicate.
        String peak_type
        File chrsz            # 2-col chromosome sizes file
    }

    command {
        set -e
        python3 $(which encode_task_reproducibility.py) \
            ${sep=' ' peaks} \
            --peaks-pr ${sep=' ' peaks_pr} \
            ${'--peak-ppr '+ peak_ppr} \
            --prefix ${prefix} \
            ${'--peak-type ' + peak_type} \
            ${'--chrsz ' + chrsz}
    }
    output {
        File optimal_peak = glob('*optimal_peak.*.gz')[0]
        File optimal_peak_bb = glob('*optimal_peak.*.bb')[0]
        File optimal_peak_starch = glob('*optimal_peak.*.starch')[0]
        File optimal_peak_hammock = glob('*optimal_peak.*.hammock.gz*')[0]
        File optimal_peak_hammock_tbi = glob('*optimal_peak.*.hammock.gz*')[1]
        File conservative_peak = glob('*conservative_peak.*.gz')[0]
        File conservative_peak_bb = glob('*conservative_peak.*.bb')[0]
        File conservative_peak_starch = glob('*conservative_peak.*.starch')[0]
        File conservative_peak_hammock = glob('*conservative_peak.*.hammock.gz*')[0]
        File conservative_peak_hammock_tbi = glob('*conservative_peak.*.hammock.gz*')[1]
        File reproducibility_qc = glob('*reproducibility.qc')[0]
        # QC metrics for optimal peak
        File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
        File peak_region_size_plot = glob('*.peak_region_size.png')[0]
        File num_peak_qc = glob('*.num_peak.qc')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 1
        disks : 'local-disk 50 SSD'
    }    
}

task gc_bias {
    input {
        File? nodup_bam
        File ref_fa

        String? picard_java_heap
    }
    Float mem_factor = 0.3
    Float input_file_size_gb = size(nodup_bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_gc_bias.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--ref-fa ' + ref_fa} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    output {
        File gc_plot = glob('*.gc_plot.png')[0]
        File gc_log = glob('*.gc.txt')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 6
        disks : 'local-disk 150 SSD'
    }
}

task qc_report {
    input {
        # optional metadata
        String pipeline_ver
        String title # name of sample
        String description # description for sample
        String? genome
        #String? encode_accession_id    # ENCODE accession ID of sample
        # workflow params
        Array[Boolean] paired_ends
        Array[Boolean] ctl_paired_ends
        String pipeline_type
        String aligner
        String peak_caller
        Int cap_num_peak
        Float idr_thresh
        Float pval_thresh
        Int xcor_trim_bp
        Int xcor_subsample_reads
        # QCs
        Array[File] samstat_qcs
        Array[File] nodup_samstat_qcs
        Array[File] dup_qcs
        Array[File] lib_complexity_qcs
        Array[File] ctl_samstat_qcs
        Array[File] ctl_nodup_samstat_qcs
        Array[File] ctl_dup_qcs
        Array[File] ctl_lib_complexity_qcs
        Array[File] xcor_plots
        Array[File] xcor_scores
        File? jsd_plot
        Array[File]? jsd_qcs
        Array[File] idr_plots
        Array[File]? idr_plots_pr
        File? idr_plot_ppr
        Array[File] frip_qcs
        Array[File] frip_qcs_pr1
        Array[File] frip_qcs_pr2
        File? frip_qc_pooled
        File? frip_qc_ppr1 
        File? frip_qc_ppr2 
        Array[File] frip_idr_qcs
        Array[File]? frip_idr_qcs_pr
        File? frip_idr_qc_ppr 
        Array[File] frip_overlap_qcs
        Array[File]? frip_overlap_qcs_pr
        File? frip_overlap_qc_ppr
        File? idr_reproducibility_qc
        File? overlap_reproducibility_qc

        Array[File] gc_plots

        Array[File] peak_region_size_qcs
        Array[File] peak_region_size_plots
        Array[File] num_peak_qcs

        File? idr_opt_peak_region_size_qc
        File? idr_opt_peak_region_size_plot
        File? idr_opt_num_peak_qc

        File? overlap_opt_peak_region_size_qc
        File? overlap_opt_peak_region_size_plot
        File? overlap_opt_num_peak_qc

        File? qc_json_ref
    }

    command {
        set -e
        python3 $(which encode_task_qc_report.py) \
            ${'--pipeline-ver ' + pipeline_ver} \
            ${"--title '" + sub(title,"'","_") + "'"} \
            ${"--desc '" + sub(description,"'","_") + "'"} \
            ${'--genome ' + genome} \
            ${'--multimapping ' + 0} \
            --paired-ends ${sep=' ' paired_ends} \
            --ctl-paired-ends ${sep=' ' ctl_paired_ends} \
            --pipeline-type ${pipeline_type} \
            --aligner ${aligner} \
            --peak-caller ${peak_caller} \
            ${'--cap-num-peak ' + cap_num_peak} \
            --idr-thresh ${idr_thresh} \
            --pval-thresh ${pval_thresh} \
            --xcor-trim-bp ${xcor_trim_bp} \
            --xcor-subsample-reads ${xcor_subsample_reads} \
            --samstat-qcs ${sep='_:_' samstat_qcs} \
            --nodup-samstat-qcs ${sep='_:_' nodup_samstat_qcs} \
            --dup-qcs ${sep='_:_' dup_qcs} \
            --lib-complexity-qcs ${sep='_:_' lib_complexity_qcs} \
            --xcor-plots ${sep='_:_' xcor_plots} \
            --xcor-scores ${sep='_:_' xcor_scores} \
            --idr-plots ${sep='_:_' idr_plots} \
            --idr-plots-pr ${sep='_:_' idr_plots_pr} \
            --ctl-samstat-qcs ${sep='_:_' ctl_samstat_qcs} \
            --ctl-nodup-samstat-qcs ${sep='_:_' ctl_nodup_samstat_qcs} \
            --ctl-dup-qcs ${sep='_:_' ctl_dup_qcs} \
            --ctl-lib-complexity-qcs ${sep='_:_' ctl_lib_complexity_qcs} \
            ${'--jsd-plot ' + jsd_plot} \
            --jsd-qcs ${sep='_:_' jsd_qcs} \
            ${'--idr-plot-ppr ' + idr_plot_ppr} \
            --frip-qcs ${sep='_:_' frip_qcs} \
            --frip-qcs-pr1 ${sep='_:_' frip_qcs_pr1} \
            --frip-qcs-pr2 ${sep='_:_' frip_qcs_pr2} \
            ${'--frip-qc-pooled ' + frip_qc_pooled} \
            ${'--frip-qc-ppr1 ' + frip_qc_ppr1} \
            ${'--frip-qc-ppr2 ' + frip_qc_ppr2} \
            --frip-idr-qcs ${sep='_:_' frip_idr_qcs} \
            --frip-idr-qcs-pr ${sep='_:_' frip_idr_qcs_pr} \
            ${'--frip-idr-qc-ppr ' + frip_idr_qc_ppr} \
            --frip-overlap-qcs ${sep='_:_' frip_overlap_qcs} \
            --frip-overlap-qcs-pr ${sep='_:_' frip_overlap_qcs_pr} \
            ${'--frip-overlap-qc-ppr ' + frip_overlap_qc_ppr} \
            ${'--idr-reproducibility-qc ' + idr_reproducibility_qc} \
            ${'--overlap-reproducibility-qc ' + overlap_reproducibility_qc} \
            --gc-plots ${sep='_:_' gc_plots} \
            --peak-region-size-qcs ${sep='_:_' peak_region_size_qcs} \
            --peak-region-size-plots ${sep='_:_' peak_region_size_plots} \
            --num-peak-qcs ${sep='_:_' num_peak_qcs} \
            ${'--idr-opt-peak-region-size-qc ' + idr_opt_peak_region_size_qc} \
            ${'--idr-opt-peak-region-size-plot ' + idr_opt_peak_region_size_plot} \
            ${'--idr-opt-num-peak-qc ' + idr_opt_num_peak_qc} \
            ${'--overlap-opt-peak-region-size-qc ' + overlap_opt_peak_region_size_qc} \
            ${'--overlap-opt-peak-region-size-plot ' + overlap_opt_peak_region_size_plot} \
            ${'--overlap-opt-num-peak-qc ' + overlap_opt_num_peak_qc} \
            --out-qc-html qc.html \
            --out-qc-json qc.json \
            ${'--qc-json-ref ' + qc_json_ref}
    }
    output {
        File report = glob('*qc.html')[0]
        File qc_json = glob('*qc.json')[0]
        Boolean qc_json_ref_match = read_string('qc_json_ref_match.txt')=='True'
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 1
        disks : 'local-disk 50 SSD'
    }    
}

### workflow system tasks
task read_genome_tsv {
    input {
        File? genome_tsv
        String? null_s
    }
    command <<<
        echo "$(basename ~{genome_tsv})" > genome_name
        # create empty files for all entries
        touch ref_fa bowtie2_idx_tar bwa_idx_tar chrsz gensz blacklist blacklist2
        touch mito_chr_name
        touch regex_bfilt_peak_chr_name

        python <<CODE
        import os
        with open('~{genome_tsv}','r') as fp:
            for line in fp:
                arr = line.strip('\n').split('\t')
                if arr:
                    key, val = arr
                    with open(key,'w') as fp2:
                        fp2.write(val)
        CODE
    >>>
    output {
        String? genome_name = read_string('genome_name')
        String? ref_fa = if size('ref_fa')==0 then null_s else read_string('ref_fa')
        String? bwa_idx_tar = if size('bwa_idx_tar')==0 then null_s else read_string('bwa_idx_tar')
        String? bowtie2_idx_tar = if size('bowtie2_idx_tar')==0 then null_s else read_string('bowtie2_idx_tar')
        String? chrsz = if size('chrsz')==0 then null_s else read_string('chrsz')
        String? gensz = if size('gensz')==0 then null_s else read_string('gensz')
        String? blacklist = if size('blacklist')==0 then null_s else read_string('blacklist')
        String? blacklist2 = if size('blacklist2')==0 then null_s else read_string('blacklist2')
        String? mito_chr_name = if size('mito_chr_name')==0 then null_s else read_string('mito_chr_name')
        String? regex_bfilt_peak_chr_name = if size('regex_bfilt_peak_chr_name')==0 then 'chr[\\dXY]+'
            else read_string('regex_bfilt_peak_chr_name')
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 1
        disks : 'local-disk 10 SSD'        
    }
}

task rounded_mean {
    input {
        Array[Int] ints
    }
    command <<<
        python <<CODE
        arr = [~{sep=',' ints}]
        with open('tmp.txt','w') as fp:
            if len(arr):
                sum_ = sum(arr)
                mean_ = sum(arr)/float(len(arr))
                fp.write('{}'.format(int(round(mean_))))
            else:
                fp.write('0')
        CODE
    >>>
    output {
        Int rounded_mean = read_int('tmp.txt')
    }
    runtime {
        cpu : 1
        memory : '2 GB'
        time : 1
        disks : 'local-disk 10 SSD'
    }
}

task raise_exception {
    input {
        String msg
        Array[String]? vals
    }
    command {
        echo -e "\n* Error: ${msg}\n" >&2
        echo -e "* Vals: ${sep=',' vals}\n" >&2
        exit 2
    }
    output {
        String error_msg = '${msg}'
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 1
        disks : 'local-disk 10 SSD'
    }
}
