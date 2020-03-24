# ENCODE TF/Histone ChIP-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

#CAPER docker quay.io/encode-dcc/chip-seq-pipeline:v1.3.7
#CAPER singularity docker://quay.io/encode-dcc/chip-seq-pipeline:v1.3.7
#CROO out_def https://storage.googleapis.com/encode-pipeline-output-definition/chip.croo.v4.json

workflow chip {
	String pipeline_ver = 'v1.3.7'
	### sample name, description
	String title = 'Untitled'
	String description = 'No description'

	# endedness for input data
	Boolean? paired_end				# to define endedness for all replciates
									#	if defined, this will override individual endedness below
	Array[Boolean] paired_ends = []	# to define endedness for individual replicate
	Boolean? ctl_paired_end
	Array[Boolean] ctl_paired_ends = []

	### mandatory genome param
	File? genome_tsv 				# reference genome data TSV file including
									# all genome-specific file paths and parameters
	# individual genome parameters
	String? genome_name				# genome name
	File? ref_fa					# reference fasta (*.fa.gz)
	File? bwa_idx_tar 				# bwa index tar (uncompressed .tar)
	File? bowtie2_idx_tar 			# bowtie2 index tar (uncompressed .tar)
	File? custom_aligner_idx_tar 	# custom aligner's index tar (uncompressed .tar)
	File? chrsz 					# 2-col chromosome sizes file
	File? blacklist 				# blacklist BED (peaks overlapping will be filtered out)
	File? blacklist2 				# 2nd blacklist (will be merged with 1st one)
	String? mito_chr_name
	String? regex_bfilt_peak_chr_name
	String? gensz 					# genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)
	File? tss 						# TSS BED file
	File? dnase 					# open chromatin region BED file
	File? prom 						# promoter region BED file
	File? enh 						# enhancer region BED file
	File? reg2map 					# file with cell type signals
	File? reg2map_bed 				# file of regions used to generate reg2map signals
	File? roadmap_meta 				# roadmap metedata

	### pipeline type
	String pipeline_type  			# tf or histone chip-seq

	String aligner = 'bowtie2'
	File? custom_align_py 			# custom align python script

	String? peak_caller 			# default: (spp for tf) and (macs2 for histone)
									# this will override the above defaults
	String? peak_type 				# default: narrowPeak for macs2, regionPeak for spp
									# this will override the above defaults
	File? custom_call_peak_py 		# custom call_peak python script

	## parameters for alignment
	Boolean align_only = false 		# disable all post-align analysis (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses involving pseudo replicates (including overlap/idr)
	Boolean enable_count_signal_track = false 		# generate count signal track
	Boolean enable_jsd = true 		# enable JSD plot generation (deeptools fingerprint)
	Boolean enable_gc_bias = true

	# parameters for aligner and filter
	Boolean use_bwa_mem_for_pe = false # THIS IS EXPERIMENTAL and BWA ONLY (use bwa mem instead of bwa aln/sam)
									# available only for PE dataset with READ_LEN>=70bp
	Int crop_length = 0 			# crop reads in FASTQs with Trimmomatic (0 by default, i.e. disabled)
	Int crop_length_tol = 2 	# keep shorter reads around crop_length
	Int xcor_trim_bp = 50 			# for cross-correlation analysis only (R1 of paired-end fastqs)
	Boolean use_filt_pe_ta_for_xcor = false # PE only. use filtered PE BAM for cross-corr.
	String dup_marker = 'picard'	# picard, sambamba
	Boolean no_dup_removal = false	# keep all dups in final BAM
	Int? mapq_thresh 				# threshold for low MAPQ reads removal
	Int mapq_thresh_bwa = 30
	Int mapq_thresh_bowtie2 = 30
	Array[String] filter_chrs = [] 	# array of chromosomes to be removed from nodup/filt BAM
									# chromosomes will be removed from both BAM header/contents
									# e.g. ['chrM', 'MT']
	Int subsample_reads = 0			# number of reads to subsample TAGALIGN
									# 0 for no subsampling. this affects all downstream analysis
	Int ctl_subsample_reads = 0		# number of reads to subsample control TAGALIGN
	Int xcor_subsample_reads = 15000000 # subsample TAG-ALIGN for xcor only (not used for other downsteam analyses)
	Int xcor_exclusion_range_min = -500
	Int? xcor_exclusion_range_max

	# parameters for peak calling
	Boolean always_use_pooled_ctl = false # always use pooled control for all exp rep.
	Float ctl_depth_ratio = 1.2 	# if ratio between controls is higher than this
									# then always use pooled control for all exp rep.
	Int? cap_num_peak
	Int cap_num_peak_spp = 300000	# cap number of raw peaks called from SPP
	Int cap_num_peak_macs2 = 500000	# cap number of raw peaks called from MACS2
	Float pval_thresh = 0.01		# p.value threshold (for MACS2 peak caller only)
	Float fdr_thresh = 0.01			# FDR threshold (for SPP peak caller only: Rscript run_spp.R -fdr)
	Float idr_thresh = 0.05			# IDR threshold

	### resources
	Int align_cpu = 4
	Int align_mem_mb = 20000
	Int align_time_hr = 48
	String align_disks = 'local-disk 400 HDD'

	Int filter_cpu = 2
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = 'local-disk 400 HDD'

	Int bam2ta_cpu = 2
	Int bam2ta_mem_mb = 10000
	Int bam2ta_time_hr = 6
	String bam2ta_disks = 'local-disk 100 HDD'

	Int spr_mem_mb = 16000

	Int jsd_cpu = 2
	Int jsd_mem_mb = 12000
	Int jsd_time_hr = 6
	String jsd_disks = 'local-disk 200 HDD'

	Int xcor_cpu = 2
	Int xcor_mem_mb = 16000	
	Int xcor_time_hr = 24
	String xcor_disks = 'local-disk 100 HDD'

	Int macs2_signal_track_mem_mb = 16000
	Int macs2_signal_track_time_hr = 24
	String macs2_signal_track_disks = 'local-disk 400 HDD'

	Int call_peak_cpu = 2
	Int call_peak_mem_mb = 16000
	Int call_peak_time_hr = 72
	String call_peak_disks = 'local-disk 200 HDD'

	String? align_trimmomatic_java_heap
	String? filter_picard_java_heap
	String? gc_bias_picard_java_heap

	#### input file definition
	# pipeline can start from any type of inputs and then leave all other types undefined
	# supported types: fastq, bam, nodup_bam (filtered bam), ta (tagAlign), peak
	# define up to 4 replicates
	# [rep_id] is for each replicate

 	### fastqs
	Array[File] fastqs_rep1_R1 = []		# FASTQs to be merged for rep1 R1
	Array[File] fastqs_rep1_R2 = [] 	# do not define if single-ended
	Array[File] fastqs_rep2_R1 = [] 	# do not define if unreplicated
	Array[File] fastqs_rep2_R2 = []		# ...
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

	Array[File] ctl_fastqs_rep1_R1 = []		# Control FASTQs to be merged for rep1 R1
	Array[File] ctl_fastqs_rep1_R2 = [] 	# do not define if single-ended
	Array[File] ctl_fastqs_rep2_R1 = [] 	# do not define if unreplicated
	Array[File] ctl_fastqs_rep2_R2 = []		# ...
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

	### other input types (bam, nodup_bam, ta)
	Array[File?] bams = [] 			# [rep_id]
	Array[File?] ctl_bams = [] 		# [rep_id]
	Array[File?] nodup_bams = [] 	# [rep_id]
	Array[File?] ctl_nodup_bams = [] # [rep_id]
	Array[File?] tas = []			# [rep_id]
	Array[File?] ctl_tas = []		# [rep_id]

	### other input types (peak)
	Array[Int?] fraglen = [] 	# [rep_id]. fragment length if inputs are peaks
	Array[File?] peaks = []		# [PAIR(rep_id1,rep_id2)]. example for 3 reps: [rep1_rep2, rep1_rep3, rep2_rep3]
	Array[File?] peaks_pr1 = []	# [rep_id]. do not define if true_rep=true
	Array[File?] peaks_pr2 = []	# [rep_id]. do not define if true_rep=true
	File? peak_ppr1				# do not define if you have a single replicate or true_rep=true
	File? peak_ppr2				# do not define if you have a single replicate or true_rep=true
	File? peak_pooled			# do not define if you have a single replicate or true_rep=true

	####################### pipeline starts here #######################
	# DO NOT DEFINE ANY VARIABLES DECLARED BELOW IN AN INPUT JSON FILE #
	# THEY ARE TEMPORARY/INTERMEDIATE SYSTEM VARIABLES                 #
	####################### pipeline starts here #######################

	# read genome data and paths
	if ( defined(genome_tsv) ) {
		call read_genome_tsv { input: genome_tsv = genome_tsv }
	}
	File? ref_fa_ = if defined(ref_fa) then ref_fa
		else read_genome_tsv.ref_fa
	File? bwa_idx_tar_ = if defined(bwa_idx_tar) then bwa_idx_tar
		else read_genome_tsv.bwa_idx_tar
	File? bowtie2_idx_tar_ = if defined(bowtie2_idx_tar) then bowtie2_idx_tar
		else read_genome_tsv.bowtie2_idx_tar
	File? custom_aligner_idx_tar_ = if defined(custom_aligner_idx_tar) then custom_aligner_idx_tar
		else read_genome_tsv.custom_aligner_idx_tar
	File? chrsz_ = if defined(chrsz) then chrsz
		else read_genome_tsv.chrsz
	String? gensz_ = if defined(gensz) then gensz
		else read_genome_tsv.gensz
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
	String? mito_chr_name_ = if defined(mito_chr_name) then mito_chr_name
		else read_genome_tsv.mito_chr_name		
	String? regex_bfilt_peak_chr_name_ = if defined(regex_bfilt_peak_chr_name) then regex_bfilt_peak_chr_name
		else read_genome_tsv.regex_bfilt_peak_chr_name
	String? genome_name_ = if defined(genome_name) then genome_name
		else if defined(read_genome_tsv.genome_name) then read_genome_tsv.genome_name
		else basename(select_first([genome_tsv, ref_fa_, chrsz_, 'None']))

	# read additional annotation data
	File? tss_ = if defined(tss) then tss
		else read_genome_tsv.tss
	File? dnase_ = if defined(dnase) then dnase
		else read_genome_tsv.dnase
	File? prom_ = if defined(prom) then prom
		else read_genome_tsv.prom
	File? enh_ = if defined(enh) then enh
		else read_genome_tsv.enh
	File? reg2map_ = if defined(reg2map) then reg2map
		else read_genome_tsv.reg2map
	File? reg2map_bed_ = if defined(reg2map_bed) then reg2map_bed
		else read_genome_tsv.reg2map_bed
	File? roadmap_meta_ = if defined(roadmap_meta) then roadmap_meta
		else read_genome_tsv.roadmap_meta

	### temp vars (do not define these)
	String aligner_ = if defined(custom_align_py) then 'custom' else aligner
	String peak_caller_ = if defined(custom_call_peak_py) then 'custom'
						else if pipeline_type=='tf' then select_first([peak_caller, 'spp'])
						else select_first([peak_caller, 'macs2'])
	String peak_type_ = if peak_caller_=='spp' then select_first([peak_type, 'regionPeak'])
						else if peak_caller_=='macs2' then select_first([peak_type, 'narrowPeak'])
						else select_first([peak_type, 'narrowPeak'])
	Boolean enable_idr = pipeline_type=='tf' # enable_idr for TF chipseq only
	String idr_rank_ = if peak_caller_=='spp' then 'signal.value'
						else if peak_caller_=='macs2' then 'p.value'
						else 'p.value'
	Int cap_num_peak_ = if peak_caller_ == 'spp' then select_first([cap_num_peak, cap_num_peak_spp])
		else select_first([cap_num_peak, cap_num_peak_macs2])
	Int mapq_thresh_ = if aligner=='bowtie2' then select_first([mapq_thresh, mapq_thresh_bowtie2])
						else select_first([mapq_thresh, mapq_thresh_bwa])

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
	# 	WDLic implementation of max(A,B,C,...)
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
	if ( !defined(chrsz_) ) {
		call raise_exception as error_genome_database { input:
			msg = 'No genome database found in your input JSON. Did you define "chip.genome_tsv" correctly?'
		}
	}
	if ( !align_only && peak_caller_ == 'spp' && num_ctl == 0 ) {
		call raise_exception as error_control_required { input:
			msg = 'SPP requires control inputs. Define control input files ("chip.ctl_*") in an input JSON file.'
		}
	}

	# align each replicate
	scatter(i in range(num_rep)) {
		# to override endedness definition for individual replicate
		# 	paired_end will override paired_ends[i]
		Boolean? paired_end_ = if !defined(paired_end) && i<length(paired_ends) then paired_ends[i]
			else paired_end

		Boolean has_input_of_align = i<length(fastqs_R1) && length(fastqs_R1[i])>0
		Boolean has_output_of_align = i<length(bams) && defined(bams[i])
		if ( has_input_of_align && !has_output_of_align ) {
			call align { input :
				fastqs_R1 = fastqs_R1[i],
				fastqs_R2 = fastqs_R2[i],
				crop_length = crop_length,
				crop_length_tol = crop_length_tol,

				aligner = aligner_,
				mito_chr_name = mito_chr_name_,
				custom_align_py = custom_align_py,
				idx_tar = if aligner=='bwa' then bwa_idx_tar_
					else if aligner=='bowtie2' then bowtie2_idx_tar_
					else custom_aligner_idx_tar_,
				paired_end = paired_end_,
				use_bwa_mem_for_pe = use_bwa_mem_for_pe,

				trimmomatic_java_heap = align_trimmomatic_java_heap,
				cpu = align_cpu,
				mem_mb = align_mem_mb,
				time_hr = align_time_hr,
				disks = align_disks,
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
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh_,
				filter_chrs = filter_chrs,
				chrsz = chrsz_,
				no_dup_removal = no_dup_removal,
				mito_chr_name = mito_chr_name_,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
				picard_java_heap = filter_picard_java_heap,
				time_hr = filter_time_hr,
				disks = filter_disks,
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
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
			}
		}
		File? ta_ = if has_output_of_bam2ta then tas[i] else bam2ta.ta

		Boolean has_input_of_spr = has_output_of_bam2ta || defined(bam2ta.ta)
		if ( has_input_of_spr && !align_only && !true_rep_only ) {
			call spr { input :
				ta = ta_,
				paired_end = paired_end_,
				mem_mb = spr_mem_mb,
			}
		}

		Boolean has_input_of_count_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
		if ( has_input_of_count_signal_track && enable_count_signal_track ) {
			# generate count signal track
			call count_signal_track { input :
				ta = ta_,
				chrsz = chrsz_,
			}
		}

		if ( enable_gc_bias && defined(nodup_bam_) && defined(ref_fa_) ) {
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

				aligner = aligner_,
				mito_chr_name = mito_chr_name_,
				custom_align_py = custom_align_py,
				idx_tar = if aligner=='bwa' then bwa_idx_tar_
					else if aligner=='bowtie2' then bowtie2_idx_tar_
					else custom_aligner_idx_tar_,
				paired_end = false,
				use_bwa_mem_for_pe = use_bwa_mem_for_pe,

				cpu = align_cpu,
				mem_mb = align_mem_mb,
				time_hr = align_time_hr,
				disks = align_disks,
			}
			# no bam deduping for xcor
			call filter as filter_R1 { input :
				bam = align_R1.bam,
				paired_end = false,
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh_,
				filter_chrs = filter_chrs,
				chrsz = chrsz_,
				no_dup_removal = true,
				mito_chr_name = mito_chr_name_,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
				picard_java_heap = filter_picard_java_heap,
				time_hr = filter_time_hr,
				disks = filter_disks,
			}
			call bam2ta as bam2ta_no_dedup_R1 { input :
				bam = filter_R1.nodup_bam,  # it's named as nodup bam but it's not deduped but just filtered
				paired_end = false,
				subsample = 0,
				mito_chr_name = mito_chr_name_,

				cpu = bam2ta_cpu,
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
			}
		}

		# special trimming/mapping for xcor (when starting from BAMs)
		Boolean has_input_of_bam2ta_no_dedup = (has_output_of_align || defined(align.bam))
			&& !defined(bam2ta_no_dedup_R1.ta)
		if ( has_input_of_bam2ta_no_dedup ) {
			call filter as filter_no_dedup { input :
				bam = bam_,
				paired_end = paired_end_,
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh_,
				filter_chrs = filter_chrs,
				chrsz = chrsz_,
				no_dup_removal = true,
				mito_chr_name = mito_chr_name_,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
				picard_java_heap = filter_picard_java_heap,
				time_hr = filter_time_hr,
				disks = filter_disks,
			}
			call bam2ta as bam2ta_no_dedup { input :
				bam = filter_no_dedup.nodup_bam,  # output name is nodup but it's not deduped
				paired_end = paired_end_,
				subsample = 0,
				mito_chr_name = mito_chr_name_,

				cpu = bam2ta_cpu,
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
			}
		}

		# use trimmed/unfilitered R1 tagAlign for paired end dataset 		
		# if not starting from fastqs, keep using old method
		#  (mapping with both ends for tag-aligns to be used for xcor)
		# subsample tagalign (non-mito) and cross-correlation analysis
		File? ta_xcor = if defined(bam2ta_no_dedup_R1.ta) then bam2ta_no_dedup_R1.ta
			else if defined(bam2ta_no_dedup.ta) then bam2ta_no_dedup.ta
			else ta_
		Boolean? paired_end_xcor = if defined(bam2ta_no_dedup_R1.ta) then false
			else paired_end_

		Boolean has_input_of_xcor = defined(ta_xcor)
		if ( has_input_of_xcor ) {
			call xcor { input :
				ta = ta_xcor,
				paired_end = paired_end_xcor,
				subsample = xcor_subsample_reads,
				mito_chr_name = mito_chr_name_,
				chip_seq_type = pipeline_type,
				exclusion_range_min = xcor_exclusion_range_min,
				exclusion_range_max = xcor_exclusion_range_max,
				cpu = xcor_cpu,
				mem_mb = xcor_mem_mb,
				time_hr = xcor_time_hr,
				disks = xcor_disks,
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
		# 	ctl_paired_end will override ctl_paired_ends[i]
		Boolean? ctl_paired_end_ = if !defined(ctl_paired_end) && i<length(ctl_paired_ends) then ctl_paired_ends[i]
			else if defined(ctl_paired_end) then ctl_paired_end
			else paired_end

		Boolean has_input_of_align_ctl = i<length(ctl_fastqs_R1) && length(ctl_fastqs_R1[i])>0
		Boolean has_output_of_align_ctl = i<length(ctl_bams) && defined(ctl_bams[i])
		if ( has_input_of_align_ctl && !has_output_of_align_ctl ) {
			call align as align_ctl { input :
				fastqs_R1 = ctl_fastqs_R1[i],
				fastqs_R2 = ctl_fastqs_R2[i],
				crop_length = crop_length,
				crop_length_tol = crop_length_tol,

				aligner = aligner_,
				mito_chr_name = mito_chr_name_,
				custom_align_py = custom_align_py,
				idx_tar = if aligner=='bwa' then bwa_idx_tar_
					else if aligner=='bowtie2' then bowtie2_idx_tar_
					else custom_aligner_idx_tar_,
				paired_end = ctl_paired_end_,
				use_bwa_mem_for_pe = use_bwa_mem_for_pe,

				trimmomatic_java_heap = align_trimmomatic_java_heap,
				cpu = align_cpu,
				mem_mb = align_mem_mb,
				time_hr = align_time_hr,
				disks = align_disks,
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
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh_,
				filter_chrs = filter_chrs,
				chrsz = chrsz_,
				no_dup_removal = no_dup_removal,
				mito_chr_name = mito_chr_name_,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
				picard_java_heap = filter_picard_java_heap,
				time_hr = filter_time_hr,
				disks = filter_disks,
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
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
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
	if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
		# pool tagaligns from pseudo replicate 1
		call pool_ta as pool_ta_pr1 { input :
			tas = spr.ta_pr1,
			prefix = 'rep-pr1',
		}
	}

	# if there are pr2 TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta_pr2 = length(select_all(spr.ta_pr2))==num_rep
	if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
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
	if ( has_input_of_count_signal_track_pooled && enable_count_signal_track && num_rep>1 ) {
		call count_signal_track as count_signal_track_pooled { input :
			ta = pool_ta.ta_pooled,
			chrsz = chrsz_,
		}
	}

	Boolean has_input_of_jsd = defined(blacklist_) && length(select_all(nodup_bam_))==num_rep
	if ( has_input_of_jsd && num_rep > 0 && enable_jsd ) {
		# fingerprint and JS-distance plot
		call jsd { input :
			nodup_bams = nodup_bam_,
			ctl_bams = ctl_nodup_bam_, # use first control only
			blacklist = blacklist_,
			mapq_thresh = mapq_thresh_,

			cpu = jsd_cpu,
			mem_mb = jsd_mem_mb,
			time_hr = jsd_time_hr,
			disks = jsd_disks,
		}
	}

	Boolean has_all_input_of_choose_ctl = length(select_all(ta_))==num_rep
		&& length(select_all(ctl_ta_))==num_ctl && num_ctl > 0
	if ( has_all_input_of_choose_ctl && !align_only ) {
		# choose appropriate control for each exp IP replicate
		# outputs:
		# 	choose_ctl.idx : control replicate index for each exp replicate 
		#					-1 means pooled ctl replicate
		call choose_ctl { input:
			tas = ta_,
			ctl_tas = ctl_ta_,
			ta_pooled = pool_ta.ta_pooled,
			ctl_ta_pooled = pool_ta_ctl.ta_pooled,
			always_use_pooled_ctl = always_use_pooled_ctl,
			ctl_depth_ratio = ctl_depth_ratio,
		}
	}

	scatter(i in range(num_rep)) {
		# make control ta array [[1,2,3,4]] -> [[1],[2],[3],[4]]
		# chosen_ctl_ta_id
		# 	>=0: control TA index (this means that control TA with this index exists)
		# 	-1: use pooled control
		#	-2: there is no control
		Int chosen_ctl_ta_id = if has_all_input_of_choose_ctl && !align_only then
			select_first([choose_ctl.chosen_ctl_ta_ids])[i] else -2
		Array[File] chosen_ctl_tas = if chosen_ctl_ta_id == -2 then []
			else if chosen_ctl_ta_id == -1 then [ select_first([pool_ta_ctl.ta_pooled]) ]
			else [ select_first([ctl_ta_[ chosen_ctl_ta_id ]]) ]
	}

	# workaround for dx error (Unsupported combination: womType: Int womValue: ([225], Array[Int]))
	Array[Int] fraglen_tmp = select_all(fraglen_)

	# we have all tas and ctl_tas (optional for histone chipseq) ready, let's call peaks
	scatter(i in range(num_rep)) {
		Boolean has_input_of_call_peak = defined(ta_[i])
		Boolean has_output_of_call_peak = i<length(peaks) && defined(peaks[i])
		if ( has_input_of_call_peak && !has_output_of_call_peak && !align_only ) {
			call call_peak { input :
				peak_caller = peak_caller_,
				peak_type = peak_type_,
				custom_call_peak_py = custom_call_peak_py,
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
				mem_mb = call_peak_mem_mb,
				disks = call_peak_disks,
				time_hr = call_peak_time_hr,
			}
		}
		File? peak_ = if has_output_of_call_peak then peaks[i]
			else call_peak.peak

		# signal track
		if ( has_input_of_call_peak && !align_only ) {
			call macs2_signal_track { input :
				tas = flatten([[ta_[i]], chosen_ctl_tas[i]]),
				gensz = gensz_,
				chrsz = chrsz_,
				pval_thresh = pval_thresh,
				fraglen = fraglen_tmp[i],

				mem_mb = macs2_signal_track_mem_mb,
				disks = macs2_signal_track_disks,
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
				custom_call_peak_py = custom_call_peak_py,
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
				mem_mb = call_peak_mem_mb,
				disks = call_peak_disks,
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
				custom_call_peak_py = custom_call_peak_py,
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
				mem_mb = call_peak_mem_mb,
				disks = call_peak_disks,
				time_hr = call_peak_time_hr,
			}
		}
		File? peak_pr2_ = if has_output_of_call_peak_pr2 then peaks_pr2[i]
			else call_peak_pr2.peak
	}

	# if ( !align_only && num_rep > 1 ) {
	# rounded mean of fragment length, which will be used for 
	#  1) calling peaks for pooled true/pseudo replicates
	#  2) calculating FRiP
	call rounded_mean as fraglen_mean { input :
		ints = fraglen_tmp,
	}
	# }

	# actually not an array
	Array[File?] chosen_ctl_ta_pooled = if !has_all_input_of_choose_ctl then []
		else if num_ctl < 2 then [ctl_ta_[0]] # choose first (only) control
		else select_all([pool_ta_ctl.ta_pooled]) # choose pooled control

	Boolean has_input_of_call_peak_pooled = defined(pool_ta.ta_pooled)
	Boolean has_output_of_call_peak_pooled = defined(peak_pooled)
	if ( has_input_of_call_peak_pooled && !has_output_of_call_peak_pooled && !align_only && num_rep>1 ) {
		# call peaks on pooled replicate
		# always call peaks for pooled replicate to get signal tracks
		call call_peak as call_peak_pooled { input :
			peak_caller = peak_caller_,
			peak_type = peak_type_,
			custom_call_peak_py = custom_call_peak_py,
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
			mem_mb = call_peak_mem_mb,
			disks = call_peak_disks,
			time_hr = call_peak_time_hr,
		}
	}
	File? peak_pooled_ = if has_output_of_call_peak_pooled then peak_pooled
		else call_peak_pooled.peak	

	# macs2 signal track for pooled rep
	if ( has_input_of_call_peak_pooled && !align_only && num_rep>1 ) {
		call macs2_signal_track as macs2_signal_track_pooled { input :
			tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz_,
			chrsz = chrsz_,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,

			mem_mb = macs2_signal_track_mem_mb,
			disks = macs2_signal_track_disks,
			time_hr = macs2_signal_track_time_hr,
		}
	}

	Boolean has_input_of_call_peak_ppr1 = defined(pool_ta_pr1.ta_pooled)
	Boolean has_output_of_call_peak_ppr1 = defined(peak_ppr1)
	if ( has_input_of_call_peak_ppr1 && !has_output_of_call_peak_ppr1 && !align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 1st pooled pseudo replicates
		call call_peak as call_peak_ppr1 { input :
			peak_caller = peak_caller_,
			peak_type = peak_type_,
			custom_call_peak_py = custom_call_peak_py,
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
			mem_mb = call_peak_mem_mb,
			disks = call_peak_disks,
			time_hr = call_peak_time_hr,
		}
	}
	File? peak_ppr1_ = if has_output_of_call_peak_ppr1 then peak_ppr1
		else call_peak_ppr1.peak

	Boolean has_input_of_call_peak_ppr2 = defined(pool_ta_pr2.ta_pooled)
	Boolean has_output_of_call_peak_ppr2 = defined(peak_ppr2)
	if ( has_input_of_call_peak_ppr2 && !has_output_of_call_peak_ppr2 && !align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 2nd pooled pseudo replicates
		call call_peak as call_peak_ppr2 { input :
			peak_caller = peak_caller_,
			peak_type = peak_type_,
			custom_call_peak_py = custom_call_peak_py,
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
			mem_mb = call_peak_mem_mb,
			disks = call_peak_disks,
			time_hr = call_peak_time_hr,
		}
	}
	File? peak_ppr2_ = if has_output_of_call_peak_ppr2 then peak_ppr2
		else call_peak_ppr2.peak

	# do IDR/overlap on all pairs of two replicates (i,j)
	# 	where i and j are zero-based indices and 0 <= i < j < num_rep
	Array[Pair[Int, Int]] pairs_ = cross(range(num_rep),range(num_rep))
	scatter( pair in pairs_ ) {
		Pair[Int, Int]? null_pair
		Pair[Int, Int]? pairs__ = if pair.left<pair.right then pair else null_pair
	}
	Array[Pair[Int, Int]] pairs = select_all(pairs__)

	if ( !align_only ) {
		scatter( pair in pairs ) {
			# pair.left = 0-based index of 1st replicate
			# pair.right = 0-based index of 2nd replicate
			# Naive overlap on every pair of true replicates
			call overlap { input :
				prefix = 'rep'+(pair.left+1)+'_vs_rep'+(pair.right+1),
				peak1 = peak_[pair.left],
				peak2 = peak_[pair.right],
				peak_pooled = peak_pooled_,
				fraglen = fraglen_mean.rounded_mean,
				peak_type = peak_type_,
				blacklist = blacklist_,
				chrsz = chrsz_,
				regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
				ta = pool_ta.ta_pooled,
			}
		}
	}

	if ( enable_idr && !align_only ) {
		scatter( pair in pairs ) {
			# pair.left = 0-based index of 1st replicate
			# pair.right = 0-based index of 2nd replicate
			# IDR on every pair of true replicates
			call idr { input :
				prefix = 'rep'+(pair.left+1)+'_vs_rep'+(pair.right+1),
				peak1 = peak_[pair.left],
				peak2 = peak_[pair.right],
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
	if ( !align_only && !true_rep_only ) {
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

	if ( !align_only && !true_rep_only && enable_idr ) {
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

	if ( !align_only && !true_rep_only && num_rep > 1 ) {
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

	if ( !align_only && !true_rep_only && num_rep > 1 ) {
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
	if ( !align_only && !true_rep_only && num_rep > 0 ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap.bfilt_overlap_peak,
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
			peak_type = peak_type_,
			chrsz = chrsz_,
		}
	}

	if ( !align_only && !true_rep_only && num_rep > 0 && enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr.bfilt_idr_peak,
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

		samstat_qcs = align.samstat_qc,
		nodup_samstat_qcs = filter.samstat_qc,
		dup_qcs = filter.dup_qc,
		lib_complexity_qcs = filter.lib_complexity_qc,
		xcor_plots = xcor.plot_png,
		xcor_scores = xcor.score,

		ctl_samstat_qcs = align_ctl.samstat_qc,
		ctl_nodup_samstat_qcs = filter_ctl.samstat_qc,
		ctl_dup_qcs = filter_ctl.dup_qc,
		ctl_lib_complexity_qcs = filter_ctl.lib_complexity_qc,

		jsd_plot = jsd.plot,
		jsd_qcs = jsd.jsd_qcs,

		frip_qcs = call_peak.frip_qc,
		frip_qcs_pr1 = call_peak_pr1.frip_qc,
		frip_qcs_pr2 = call_peak_pr2.frip_qc,
		frip_qc_pooled = call_peak_pooled.frip_qc,
		frip_qc_ppr1 = call_peak_ppr1.frip_qc,
		frip_qc_ppr2 = call_peak_ppr2.frip_qc,

		idr_plots = idr.idr_plot,
		idr_plots_pr = idr_pr.idr_plot,
		idr_plot_ppr = idr_ppr.idr_plot,
		frip_idr_qcs = idr.frip_qc,
		frip_idr_qcs_pr = idr_pr.frip_qc,
		frip_idr_qc_ppr = idr_ppr.frip_qc,
		frip_overlap_qcs = overlap.frip_qc,
		frip_overlap_qcs_pr = overlap_pr.frip_qc,
		frip_overlap_qc_ppr = overlap_ppr.frip_qc,
		idr_reproducibility_qc = reproducibility_idr.reproducibility_qc,
		overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,

		gc_plots = gc_bias.gc_plot,

		peak_region_size_qcs = call_peak.peak_region_size_qc,
		peak_region_size_plots = call_peak.peak_region_size_plot,
		num_peak_qcs = call_peak.num_peak_qc,

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
	Array[File] fastqs_R1 		# [merge_id]
	Array[File] fastqs_R2
	Int? trim_bp			# this is for R1 only
	Int crop_length
	Int crop_length_tol
	String aligner
	String mito_chr_name
	Int? multimapping
	File? custom_align_py	
	File? idx_tar			# reference index tar
	Boolean paired_end
	Boolean use_bwa_mem_for_pe

	Float trimmomatic_java_heap_factor = 0.9
	String? trimmomatic_java_heap
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

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
				--out-dir-R1 R1$NEW_SUFFIX \
				${if paired_end then '--out-dir-R2 R2$NEW_SUFFIX' else ''} \
				${'--trimmomatic-java-heap ' + if defined(trimmomatic_java_heap) then trimmomatic_java_heap else (mem_mb * trimmomatic_java_heap_factor + 'M')} \
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
				${'--nth ' + cpu}

		elif [ '${aligner}' == 'bowtie2' ]; then
		 	python3 $(which encode_task_bowtie2.py) \
				${idx_tar} \
				R1$SUFFIX/*.fastq.gz \
				${if paired_end then 'R2$SUFFIX/*.fastq.gz' else ''} \
				${'--multimapping ' + multimapping} \
				${if paired_end then '--paired-end' else ''} \
				${'--nth ' + cpu}
		else
			python3 ${custom_align_py} \
				${idx_tar} \
				R1$SUFFIX/*.fastq.gz \
				${if paired_end then 'R2$SUFFIX/*.fastq.gz' else ''} \
				${if paired_end then '--paired-end' else ''} \
				${'--nth ' + cpu}
		fi 

		python3 $(which encode_task_post_align.py) \
			R1$SUFFIX/*.fastq.gz $(ls *.bam) \
			${'--mito-chr-name ' + mito_chr_name} \
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
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
		preemptible: 0
	}
}

task filter {
	File bam
	Boolean paired_end
	String dup_marker 			# picard.jar MarkDuplicates (picard) or 
								# sambamba markdup (sambamba)
	Int mapq_thresh				# threshold for low MAPQ reads removal
	Array[String] filter_chrs 	# chrs to be removed from final (nodup/filt) BAM
	File chrsz					# 2-col chromosome sizes file
	Boolean no_dup_removal 		# no dupe reads removal when filtering BAM
	String mito_chr_name

	Int cpu
	Int mem_mb
	Float picard_java_heap_factor = 0.9
	String? picard_java_heap
	Int time_hr
	String disks

	command {
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
			${'--nth ' + cpu} \
			${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (mem_mb * picard_java_heap_factor + 'M')}
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
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
	}
}

task bam2ta {
	File bam
	Boolean paired_end
	String mito_chr_name 		# mito chromosome name
	Int subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python3 $(which encode_task_bam2ta.py) \
			${bam} \
			--disable-tn5-shift \
			${if paired_end then '--paired-end' else ''} \
			${'--mito-chr-name ' + mito_chr_name} \
			${'--subsample ' + subsample} \
			${'--nth ' + cpu}
	}
	output {
		File ta = glob('*.tagAlign.gz')[0]
	}
	runtime {
		cpu : cpu
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
	}
}

task spr { # make two self pseudo replicates
	File ta
	Boolean paired_end

	Int mem_mb

	command {
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
		memory : '${mem_mb} MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}
}

task pool_ta {
	Array[File?] tas
	Int? col 			# number of columns in pooled TA
	String? prefix 		# basename prefix

	command {
		python3 $(which encode_task_pool_ta.py) \
			${sep=' ' tas} \
			${'--prefix ' + prefix} \
			${'--col ' + col}
	}
	output {
		File ta_pooled = glob('*.tagAlign.gz')[0]
	}
	runtime {
		cpu : 1
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}
}

task xcor {
	File ta
	Boolean paired_end
	String mito_chr_name
	Int subsample  # number of reads to subsample TAGALIGN
					# this will be used for xcor only
					# will not affect any downstream analysis
	String? chip_seq_type
	Int? exclusion_range_min
	Int? exclusion_range_max

	Int cpu
	Int mem_mb	
	Int time_hr
	String disks

	command {
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
		Int fraglen = read_int(glob('*.cc.fraglen.txt')[0])
	}
	runtime {
		cpu : cpu
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
	}
}

task jsd {
	Array[File?] nodup_bams
	Array[File?] ctl_bams
	File blacklist
	Int mapq_thresh

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python3 $(which encode_task_jsd.py) \
			${sep=' ' nodup_bams} \
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
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
	}
}

task choose_ctl {
	Array[File?] tas
	Array[File?] ctl_tas
	File? ta_pooled
	File? ctl_ta_pooled
	Boolean always_use_pooled_ctl # always use pooled control for all exp rep.
	Float ctl_depth_ratio 		# if ratio between controls is higher than this
								# then always use pooled control for all exp rep.
	command {
		python3 $(which encode_task_choose_ctl.py) \
			--tas ${sep=' ' tas} \
			--ctl-tas ${sep=' ' ctl_tas} \
			${'--ta-pooled ' + ta_pooled} \
			${'--ctl-ta-pooled ' + ctl_ta_pooled} \
			${if always_use_pooled_ctl then '--always-use-pooled-ctl' else ''} \
			${'--ctl-depth-ratio ' + ctl_depth_ratio}
	}
	output {
		Array[Int] chosen_ctl_ta_ids = read_lines('chosen_ctl.tsv')
	}
	runtime {
		cpu : 1
		memory : '8000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}	
}

task count_signal_track {
	File ta 			# tag-align
	File chrsz			# 2-col chromosome sizes file

	command {
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
		memory : '8000 MB'
		time : 4
		disks : 'local-disk 50 HDD'
	}
}

task call_peak {
	String peak_caller
	String peak_type
	File? custom_call_peak_py

	Array[File?] tas	# [ta, control_ta]. control_ta is optional
	Int fraglen 		# fragment length from xcor
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int cap_num_peak	# cap number of raw peaks called from MACS2
	Float pval_thresh 	# p.value threshold for MACS2
	Float? fdr_thresh 	# FDR threshold for SPP

	File? blacklist 	# blacklist BED to filter raw peaks
	String? regex_bfilt_peak_chr_name

	Int cpu	
	Int mem_mb
	Int time_hr
	String disks

	command {
		set -e

		if [ '${peak_caller}' == 'macs2' ]; then
			python3 $(which encode_task_macs2_chip.py) \
				${sep=' ' tas} \
				${'--gensz '+ gensz} \
				${'--chrsz ' + chrsz} \
				${'--fraglen ' + fraglen} \
				${'--cap-num-peak ' + cap_num_peak} \
				${'--pval-thresh '+ pval_thresh}

		elif [ '${peak_caller}' == 'spp' ]; then
			python3 $(which encode_task_spp.py) \
				${sep=' ' tas} \
				${'--fraglen ' + fraglen} \
				${'--cap-num-peak ' + cap_num_peak} \
				${'--fdr-thresh '+ fdr_thresh} \
				${'--nth ' + cpu}

		else
			python3 ${custom_call_peak_py} \
				${sep=' ' tas} \
				${'--gensz '+ gensz} \
				${'--chrsz ' + chrsz} \
				${'--fraglen ' + fraglen} \
				${'--cap-num-peak ' + cap_num_peak} \
				${'--pval-thresh '+ pval_thresh}
				${'--fdr-thresh '+ fdr_thresh}
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
		# generated by custom_call_peak_py
		File peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
		# generated by post_call_peak py
		File bfilt_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
		File bfilt_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
		File bfilt_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
		File bfilt_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
		File frip_qc = glob('*.frip.qc')[0]
		File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
		File peak_region_size_plot = glob('*.peak_region_size.png')[0]
		File num_peak_qc = glob('*.num_peak.qc')[0]
	}
	runtime {
		cpu : if peak_caller == 'macs2' then 1 else cpu
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
		preemptible: 0		
	}
}

task macs2_signal_track {
	Array[File?] tas	# [ta, control_ta]. control_ta is optional
	Int fraglen 		# fragment length from xcor
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Float pval_thresh 	# p.value threshold

	Int mem_mb
	Int time_hr
	String disks

	command {
		python3 $(which encode_task_macs2_signal_track_chip.py) \
			${sep=' ' tas} \
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
		memory : '${mem_mb} MB'
		time : time_hr
		disks : disks
		preemptible: 0
	}
}

task idr {
	String prefix 		# prefix for IDR output file
	File peak1 			
	File peak2
	File peak_pooled
	Float idr_thresh
	File? blacklist 	# blacklist BED to filter raw peaks
	String regex_bfilt_peak_chr_name
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
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
		File bfilt_idr_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
		File bfilt_idr_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
		File idr_plot = glob('*.txt.png')[0]
		File idr_unthresholded_peak = glob('*.txt.gz')[0]
		File idr_log = glob('*.idr*.log')[0]
		File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
	}
	runtime {
		cpu : 1
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}	
}

task overlap {
	String prefix 	# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File? blacklist # blacklist BED to filter raw peaks
	String regex_bfilt_peak_chr_name
	# parameters to compute FRiP
	File? ta		# to calculate FRiP
	Int fraglen 	# fragment length from xcor (for FRIP)
	File chrsz		# 2-col chromosome sizes file
	String peak_type

	command {
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
		File bfilt_overlap_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
		File bfilt_overlap_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
		File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
	}
	runtime {
		cpu : 1
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}	
}

task reproducibility {
	String prefix
	Array[File]? peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File?] peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.
	String peak_type
	File chrsz			# 2-col chromosome sizes file

	command {
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
		File optimal_peak_hammock = glob('*optimal_peak.*.hammock.gz*')[0]
		File optimal_peak_hammock_tbi = glob('*optimal_peak.*.hammock.gz*')[1]
		File conservative_peak = glob('*conservative_peak.*.gz')[0]
		File conservative_peak_bb = glob('*conservative_peak.*.bb')[0]
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
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}	
}

task gc_bias {
	File nodup_bam
	File ref_fa

	Int mem_mb = 10000
	Float picard_java_heap_factor = 0.9
	String? picard_java_heap

	command {
		python3 $(which encode_task_gc_bias.py) \
			${'--nodup-bam ' + nodup_bam} \
			${'--ref-fa ' + ref_fa} \
			${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (mem_mb * picard_java_heap_factor + 'M')}
	}
	output {
		File gc_plot = glob('*.gc_plot.png')[0]
		File gc_log = glob('*.gc.txt')[0]
	}
	runtime {
		cpu : 1
		memory : '${mem_mb} MB'
		time : 6
		disks : 'local-disk 100 HDD'
	}
}

# gather all outputs and generate 
# - qc.html		: organized final HTML report
# - qc.json		: all QCs
task qc_report {
	# optional metadata
	String pipeline_ver
 	String title # name of sample
	String description # description for sample
	String? genome
	#String? encode_accession_id	# ENCODE accession ID of sample
	# workflow params
	Array[Boolean?] paired_ends
	Array[Boolean?] ctl_paired_ends
	String pipeline_type
	String aligner
	String peak_caller
	Int cap_num_peak
	Float idr_thresh
	Float pval_thresh
	Int xcor_trim_bp
	Int xcor_subsample_reads
	# QCs
	Array[File?] samstat_qcs
	Array[File?] nodup_samstat_qcs
	Array[File?] dup_qcs
	Array[File?] lib_complexity_qcs
	Array[File?] ctl_samstat_qcs
	Array[File?] ctl_nodup_samstat_qcs
	Array[File?] ctl_dup_qcs
	Array[File?] ctl_lib_complexity_qcs
	Array[File?] xcor_plots
	Array[File?] xcor_scores
	File? jsd_plot
	Array[File]? jsd_qcs
	Array[File]? idr_plots
	Array[File]? idr_plots_pr
	File? idr_plot_ppr
	Array[File?] frip_qcs
	Array[File?] frip_qcs_pr1
	Array[File?] frip_qcs_pr2
	File? frip_qc_pooled
	File? frip_qc_ppr1 
	File? frip_qc_ppr2 
	Array[File]? frip_idr_qcs
	Array[File]? frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File]? frip_overlap_qcs
	Array[File]? frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	Array[File?] gc_plots

	Array[File?] peak_region_size_qcs
	Array[File?] peak_region_size_plots
	Array[File?] num_peak_qcs

	File? idr_opt_peak_region_size_qc
	File? idr_opt_peak_region_size_plot
	File? idr_opt_num_peak_qc

	File? overlap_opt_peak_region_size_qc
	File? overlap_opt_peak_region_size_plot
	File? overlap_opt_num_peak_qc

	File? qc_json_ref

	command {
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
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}	
}

### workflow system tasks
task read_genome_tsv {
	File genome_tsv

	String? null_s
	command <<<
		# create empty files for all entries
		touch genome_name
		touch ref_fa bowtie2_idx_tar bwa_idx_tar chrsz gensz blacklist blacklist2
		touch custom_aligner_idx_tar
		touch tss tss_enrich # for backward compatibility
		touch dnase prom enh reg2map reg2map_bed roadmap_meta
		touch mito_chr_name
		touch regex_bfilt_peak_chr_name

		python <<CODE
		import os
		with open('${genome_tsv}','r') as fp:
			for line in fp:
				arr = line.strip('\n').split('\t')
				if arr:
					key, val = arr
					with open(key,'w') as fp2:
						fp2.write(val)
		CODE
	>>>
	output {
		String? genome_name = if size('genome_name')==0 then basename(genome_tsv) else read_string('genome_name')
		String? ref_fa = if size('ref_fa')==0 then null_s else read_string('ref_fa')
		String? bwa_idx_tar = if size('bwa_idx_tar')==0 then null_s else read_string('bwa_idx_tar')
		String? bowtie2_idx_tar = if size('bowtie2_idx_tar')==0 then null_s else read_string('bowtie2_idx_tar')
		String? custom_aligner_idx_tar = if size('custom_aligner_idx_tar')==0 then null_s else read_string('custom_aligner_idx_tar')
		String? chrsz = if size('chrsz')==0 then null_s else read_string('chrsz')
		String? gensz = if size('gensz')==0 then null_s else read_string('gensz')
		String? blacklist = if size('blacklist')==0 then null_s else read_string('blacklist')
		String? blacklist2 = if size('blacklist2')==0 then null_s else read_string('blacklist2')
		String? mito_chr_name = if size('mito_chr_name')==0 then null_s else read_string('mito_chr_name')
		String? regex_bfilt_peak_chr_name = if size('regex_bfilt_peak_chr_name')==0 then 'chr[\\dXY]+'
			else read_string('regex_bfilt_peak_chr_name')
		# optional data
		String? tss = if size('tss')!=0 then read_string('tss')
			else if size('tss_enrich')!=0 then read_string('tss_enrich') else null_s
		String? dnase = if size('dnase')==0 then null_s else read_string('dnase')
		String? prom = if size('prom')==0 then null_s else read_string('prom')
		String? enh = if size('enh')==0 then null_s else read_string('enh')
		String? reg2map = if size('reg2map')==0 then null_s else read_string('reg2map')
		String? reg2map_bed = if size('reg2map_bed')==0 then null_s else read_string('reg2map_bed')
		String? roadmap_meta = if size('roadmap_meta')==0 then null_s else read_string('roadmap_meta')
	}
	runtime {
		maxRetries : 0
		cpu : 1
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'		
	}
}

task rounded_mean {
	Array[Int] ints
	command <<<
		python <<CODE
		arr = [${sep=',' ints}]
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
		memory : '4000 MB'
		time : 1
		disks : 'local-disk 50 HDD'
	}	
}

task raise_exception {
	String msg
	command {
		echo -e "\n* Error: ${msg}\n" >&2
		exit 2
	}
	output {
		String error_msg = '${msg}'
	}
	runtime {
		maxRetries : 0
	}
}
