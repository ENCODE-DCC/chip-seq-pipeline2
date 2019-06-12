# ENCODE DCC TF/Histone ChIP-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

#CAPER docker quay.io/encode-dcc/chip-seq-pipeline:v1.2.1
#CAPER singularity docker://quay.io/encode-dcc/chip-seq-pipeline:v1.2.1
#CROO out_def https://storage.googleapis.com/encode-pipeline-output-definition/chip.out_def.json

workflow chip {
	String pipeline_ver = 'v1.2.1'
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
	File? ref_fa					# reference fasta (*.fa.gz)
	File? bwa_idx_tar 				# bwa index tar (uncompressed)
	File? chrsz 					# 2-col chromosome sizes file
	File? blacklist 				# blacklist BED (peaks overlapping will be filtered out)
	String? gensz 					# genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)

	### pipeline type
	String pipeline_type  		# tf or histone chip-eq
	String? peak_caller 		# default: (spp for tf) and (macs2 for histone)

	### optional but important
	Boolean align_only = false 		# disable all post-align analysis (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses for pseudo replicates
									# overlap and idr will also be disabled
	Boolean disable_fingerprint = false # no JSD plot generation (deeptools fingerprint)
	Boolean enable_count_signal_track = false 		# generate count signal track
	Boolean use_bwa_mem_for_pe = false # THIS IS EXPERIMENTAL
									# use bwa mem instead of bwa aln + bwa sampe
									#for PE dataset with read len>=70bp

	Int xcor_pe_trim_bp = 50 		# for cross-correlation analysis only (R1 of paired-end fastqs)
	Boolean use_filt_pe_ta_for_xcor = false # PE only. use filtered PE BAM for cross-corr.

	String dup_marker = 'picard'	# picard.jar MarkDuplicates (picard) or 
									# sambamba markdup (sambamba)
	Int mapq_thresh = 30			# threshold for low MAPQ reads removal
	Boolean no_dup_removal = false	# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
									# filtered bam with dupes

	String mito_chr_name = 'chrM' 	# name of mito chromosome. THIS IS NOT A REG-EX! you can define only one chromosome name for mito.
	String regex_filter_reads = 'chrM' 	# Perl-style regular expression pattern for chr name to filter out reads
                        			# those reads will be excluded from peak calling
	Int subsample_reads = 0			# number of reads to subsample TAGALIGN
									# 0 for no subsampling. this affects all downstream analysis
	Int ctl_subsample_reads = 0		# number of reads to subsample control TAGALIGN

	Int xcor_subsample_reads = 15000000	# number of reads to subsample TAGALIGN
									# this will be used for xcor only
									# will not affect any downstream analysis
	Int xcor_exclusion_range_min = -500
	Int? xcor_exclusion_range_max

	Boolean keep_irregular_chr_in_bfilt_peak = false # when filtering with blacklist
								# do not filter peaks with irregular chr name
								# and just keep them in bfilt_peak file
								# (e.g. keep chr1_AABBCC, AABR07024382.1, ...)
								# reg-ex pattern for regular chr names: /chr[\dXY]+[ \t]/
	Boolean always_use_pooled_ctl = false # always use pooled control for all exp rep.
	Float ctl_depth_ratio = 1.2 	# if ratio between controls is higher than this
									# then always use pooled control for all exp rep.

	### task-specific variables but defined in workflow level (limit of WDL)
	Int macs2_cap_num_peak = 500000	# cap number of raw peaks called from MACS2
	Float pval_thresh = 0.01		# p.value threshold
	Float idr_thresh = 0.05			# IDR threshold
	Int spp_cap_num_peak = 300000	# cap number of raw peaks called from SPP

	### resources
	Int bwa_cpu = 4
	Int bwa_mem_mb = 20000
	Int bwa_time_hr = 48
	String bwa_disks = "local-disk 200 HDD"

	Int filter_cpu = 2
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = "local-disk 400 HDD"

	Int bam2ta_cpu = 2
	Int bam2ta_mem_mb = 10000
	Int bam2ta_time_hr = 6
	String bam2ta_disks = "local-disk 100 HDD"

	Int spr_mem_mb = 16000

	Int fingerprint_cpu = 2
	Int fingerprint_mem_mb = 12000
	Int fingerprint_time_hr = 6
	String fingerprint_disks = "local-disk 200 HDD"

	Int xcor_cpu = 2
	Int xcor_mem_mb = 16000	
	Int xcor_time_hr = 24
	String xcor_disks = "local-disk 100 HDD"

	Int macs2_mem_mb = 16000
	Int macs2_time_hr = 24
	String macs2_disks = "local-disk 200 HDD"

	Int spp_cpu = 2
	Int spp_mem_mb = 16000
	Int spp_time_hr = 72
	String spp_disks = "local-disk 200 HDD"

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
	Array[File?] merged_fastqs_R1 = []
	Array[File?] merged_fastqs_R2 = [] 	
	Array[File?] ctl_merged_fastqs_R1 = []
	Array[File?] ctl_merged_fastqs_R2 = [] 	
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
	File? chrsz_ = if defined(chrsz) then chrsz
		else read_genome_tsv.chrsz
	String? gensz_ = if defined(gensz) then gensz
		else read_genome_tsv.gensz
	File? blacklist_ = if defined(blacklist) then blacklist
		else read_genome_tsv.blacklist

	### temp vars (do not define these)
	String peak_caller_ = if pipeline_type=='tf' then select_first([peak_caller, 'spp'])
						else select_first([peak_caller, 'macs2'])
	String peak_type = if peak_caller_=='spp' then 'regionPeak'
						else if peak_caller_=='macs2' then 'narrowPeak'
						else 'narrowPeak'
	Boolean enable_idr = pipeline_type=='tf' # enable_idr for TF chipseq only
	String idr_rank = if peak_caller_=='spp' then 'signal.value'
						else if peak_caller_=='macs2' then 'p.value'
						else 'p.value'

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
	Int num_rep_merged_fastq = if length(merged_fastqs_R1)<num_rep_fastq then num_rep_fastq
		else length(merged_fastqs_R1)
	Int num_rep_bam = if length(bams)<num_rep_merged_fastq then num_rep_merged_fastq
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
	Int num_ctl_merged_fastq = if length(ctl_merged_fastqs_R1)<num_ctl_fastq then num_ctl_fastq
		else length(ctl_merged_fastqs_R1)
	Int num_ctl_bam = if length(ctl_bams)<num_ctl_merged_fastq then num_ctl_merged_fastq
		else length(ctl_bams)
	Int num_ctl_nodup_bam = if length(ctl_nodup_bams)<num_ctl_bam then num_ctl_bam
		else length(ctl_nodup_bams)
	Int num_ctl_ta = if length(ctl_tas)<num_ctl_nodup_bam then num_ctl_nodup_bam
		else length(ctl_tas)
	Int num_ctl = num_ctl_ta

	# align each replicate
	scatter(i in range(num_rep)) {
		# to override endedness definition for individual replicate
		# 	paired_end will override paired_ends[i]
		Boolean? paired_end_ = if !defined(paired_end) && i<length(paired_ends) then paired_ends[i]
			else paired_end
		Boolean has_input_of_merge_fastq = i<length(fastqs_R1) && length(fastqs_R1[i])>0
		Boolean has_output_of_merge_fastq = i<length(merged_fastqs_R1) &&
			defined(merged_fastqs_R1[i])
		if ( has_input_of_merge_fastq && !has_output_of_merge_fastq ) {
			# merge fastqs
			call merge_fastq { input :
				fastqs_R1 = fastqs_R1[i],
				fastqs_R2 = fastqs_R2[i],
				paired_end = paired_end_,
			}
		}
		File? merged_fastq_R1_ = if has_output_of_merge_fastq then merged_fastqs_R1[i]
			else merge_fastq.merged_fastq_R1
		File? merged_fastq_R2_ = if i<length(merged_fastqs_R2) &&
			defined(merged_fastqs_R2[i]) then merged_fastqs_R2[i]
			else merge_fastq.merged_fastq_R2

		Boolean has_input_of_bwa = has_output_of_merge_fastq || defined(merge_fastq.merged_fastq_R1)
		Boolean has_output_of_bwa = i<length(bams) && defined(bams[i])
		if ( has_input_of_bwa && !has_output_of_bwa ) {
			call bwa { input :
				bwa_idx_tar = bwa_idx_tar_,
				fastq_R1 = merged_fastq_R1_,
				fastq_R2 = merged_fastq_R2_,
				paired_end = paired_end_,
				use_bwa_mem_for_pe = use_bwa_mem_for_pe,
				cpu = bwa_cpu,
				mem_mb = bwa_mem_mb,
				time_hr = bwa_time_hr,
				disks = bwa_disks,
			}
		}
		File? bam_ = if has_output_of_bwa then bams[i] else bwa.bam

		Boolean has_input_of_filter = has_output_of_bwa || defined(bwa.bam)
		Boolean has_output_of_filter = i<length(nodup_bams) && defined(nodup_bams[i])
		# skip if we already have output of this step
		if ( has_input_of_filter && !has_output_of_filter ) {
			call filter { input :
				bam = bam_,
				paired_end = paired_end_,
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh,
				no_dup_removal = no_dup_removal,
				mito_chr_name = mito_chr_name,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
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
				regex_grep_v_ta = regex_filter_reads,
				subsample = subsample_reads,
				paired_end = paired_end_,
				mito_chr_name = mito_chr_name,

				cpu = bam2ta_cpu,
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
			}
		}
		File? ta_ = if has_output_of_bam2ta then tas[i] else bam2ta.ta

		# convert unfiltered BAM to a special TAG-ALIGN for xcor 
		Boolean has_input_of_bam2ta_no_filt = has_output_of_bwa || defined(bwa.bam)
		if ( has_input_of_bam2ta_no_filt ) {
			call bam2ta as bam2ta_no_filt { input :
				bam = bam_,
				paired_end = paired_end_,
				subsample = 0,
				regex_grep_v_ta = regex_filter_reads,
				mito_chr_name = mito_chr_name,

				cpu = bam2ta_cpu,
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
			}
		}

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

		Boolean has_input_of_trim_fastq = has_output_of_merge_fastq || defined(merge_fastq.merged_fastq_R1)
		if ( has_input_of_trim_fastq && paired_end_ ) {
			# special trimming for paired end samples (for cross-corr analysis)
			call trim_fastq { input :
				fastq = merged_fastq_R1_,
				trim_bp = xcor_pe_trim_bp,
			}
			call bwa as bwa_R1 { input :
				bwa_idx_tar = bwa_idx_tar_,
				fastq_R1 = trim_fastq.trimmed_fastq,
				paired_end = false,
				use_bwa_mem_for_pe = use_bwa_mem_for_pe,
				cpu = bwa_cpu,
				mem_mb = bwa_mem_mb,
				time_hr = bwa_time_hr,
				disks = bwa_disks,
			}
			# no bam filtering for xcor
			call bam2ta as bam2ta_no_filt_R1 { input :
				bam = bwa_R1.bam,
				paired_end = false,
				subsample = 0,
				regex_grep_v_ta = regex_filter_reads,
				mito_chr_name = mito_chr_name,

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
		File? ta_xcor = if defined(bam2ta_no_filt_R1.ta) then bam2ta_no_filt_R1.ta
			else if defined(bam2ta_no_filt.ta) then bam2ta_no_filt.ta
			else ta_
		Boolean? paired_end_xcor = if defined(bam2ta_no_filt_R1.ta) then false
			else paired_end_

		Boolean has_input_of_xcor = defined(ta_xcor)
		if ( has_input_of_xcor ) {
			call xcor { input :
				ta = ta_xcor,
				paired_end = paired_end_xcor,
				subsample = xcor_subsample_reads,
				mito_chr_name = mito_chr_name,
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

		Boolean has_input_of_merge_fastq_ctl = i<length(ctl_fastqs_R1) && length(ctl_fastqs_R1[i])>0
		Boolean has_output_of_merge_fastq_ctl = i<length(ctl_merged_fastqs_R1) &&
			defined(ctl_merged_fastqs_R1[i])
		if ( has_input_of_merge_fastq_ctl && !has_output_of_merge_fastq_ctl ) {
			# merge fastqs
			call merge_fastq as merge_fastq_ctl { input :
				fastqs_R1 = ctl_fastqs_R1[i],
				fastqs_R2 = ctl_fastqs_R2[i],
				paired_end = ctl_paired_end_,
			}
		}
		File? ctl_merged_fastq_R1_ = if has_output_of_merge_fastq_ctl then ctl_merged_fastqs_R1[i]
			else merge_fastq_ctl.merged_fastq_R1
		File? ctl_merged_fastq_R2_ = if i<length(ctl_merged_fastqs_R2) &&
			defined(ctl_merged_fastqs_R2[i]) then ctl_merged_fastqs_R2[i]
			else merge_fastq_ctl.merged_fastq_R2

		Boolean has_input_of_bwa_ctl = has_output_of_merge_fastq_ctl || defined(merge_fastq_ctl.merged_fastq_R1)
		Boolean has_output_of_bwa_ctl = i<length(ctl_bams) && defined(ctl_bams[i])
		if ( has_input_of_bwa_ctl && !has_output_of_bwa_ctl ) {
			call bwa as bwa_ctl { input :
				bwa_idx_tar = bwa_idx_tar_,
				fastq_R1 = ctl_merged_fastq_R1_,
				fastq_R2 = ctl_merged_fastq_R2_,
				paired_end = ctl_paired_end_,
				use_bwa_mem_for_pe = use_bwa_mem_for_pe,
				cpu = bwa_cpu,
				mem_mb = bwa_mem_mb,
				time_hr = bwa_time_hr,
				disks = bwa_disks,
			}
		}
		File? ctl_bam_ = if has_output_of_bwa_ctl then ctl_bams[i] else bwa_ctl.bam

		Boolean has_input_of_filter_ctl = has_output_of_bwa_ctl || defined(bwa_ctl.bam)
		Boolean has_output_of_filter_ctl = i<length(ctl_nodup_bams) && defined(ctl_nodup_bams[i])
		# skip if we already have output of this step
		if ( has_input_of_filter_ctl && !has_output_of_filter_ctl ) {
			call filter as filter_ctl { input :
				bam = ctl_bam_,
				paired_end = ctl_paired_end_,
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh,
				no_dup_removal = no_dup_removal,
				mito_chr_name = mito_chr_name,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
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
				regex_grep_v_ta = regex_filter_reads,
				subsample = subsample_reads,
				paired_end = ctl_paired_end_,
				mito_chr_name = mito_chr_name,

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
		}
	}

	# if there are pr1 TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta_pr1 = length(select_all(spr.ta_pr1))==num_rep
	if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
		# pool tagaligns from pseudo replicate 1
		call pool_ta as pool_ta_pr1 { input :
			tas = spr.ta_pr1,
		}
	}

	# if there are pr2 TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta_pr2 = length(select_all(spr.ta_pr2))==num_rep
	if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
		# pool tagaligns from pseudo replicate 2
		call pool_ta as pool_ta_pr2 { input :
			tas = spr.ta_pr2,
		}
	}

	# if there are CTL TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta_ctl = length(select_all(ctl_ta_))==num_ctl
	if ( has_all_inputs_of_pool_ta_ctl && num_ctl>1 ) {
		# pool tagaligns from true replicates
		call pool_ta as pool_ta_ctl { input :
			tas = ctl_ta_,
		}
	}

	Boolean has_input_of_count_signal_track_pooled = defined(pool_ta.ta_pooled)
	if ( has_input_of_count_signal_track_pooled && enable_count_signal_track && num_rep>1 ) {
		call count_signal_track as count_signal_track_pooled { input :
			ta = pool_ta.ta_pooled,
			chrsz = chrsz_,
		}
	}

	Boolean has_input_of_fingerprint = defined(blacklist_) && #basename(blacklist_) != 'null' &&
		length(select_all(nodup_bam_))==num_rep &&
		num_ctl>0 && defined(ctl_nodup_bam_[0])
	if ( has_input_of_fingerprint && !disable_fingerprint ) {
		# fingerprint and JS-distance plot
		call fingerprint { input :
			nodup_bams = nodup_bam_,
			ctl_bam = ctl_nodup_bam_[0], # use first control only
			blacklist = blacklist_,

			cpu = fingerprint_cpu,
			mem_mb = fingerprint_mem_mb,
			time_hr = fingerprint_time_hr,
			disks = fingerprint_disks,
		}
	}

	Boolean has_all_input_of_choose_ctl = length(select_all(ta_))==num_rep
		&& length(select_all(ctl_ta_))==num_ctl && num_ctl > 0
	if ( has_all_input_of_choose_ctl ) {
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

	# make control ta array [[1,2,3,4]] -> [[1],[2],[3],[4]], will be zipped with exp ta array latter
	Array[Array[File]] chosen_ctl_tas =
		if has_all_input_of_choose_ctl then transpose(select_all([choose_ctl.chosen_ctl_tas]))
		else [[],[],[],[],[],[],[],[],[],[]]

	# workaround for dx error (Unsupported combination: womType: Int womValue: ([225], Array[Int]))
	Array[Int] fraglen_tmp = select_all(fraglen_)

	# we have all tas and ctl_tas (optional for histone chipseq) ready, let's call peaks
	scatter(i in range(num_rep)) {
		Boolean has_input_of_peak_call = defined(ta_[i])
		Boolean has_output_of_peak_call = i<length(peaks) && defined(peaks[i])
		if ( has_input_of_peak_call && !has_output_of_peak_call &&
			 peak_caller_ == 'spp' && !align_only ) {
			call spp { input :
				tas = flatten([[ta_[i]], chosen_ctl_tas[i]]),
				chrsz = chrsz_,
				cap_num_peak = spp_cap_num_peak,
				fraglen = fraglen_tmp[i],
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
	
				cpu = spp_cpu,
				mem_mb = spp_mem_mb,
				disks = spp_disks,
				time_hr = spp_time_hr,
			}
		}
		if ( has_input_of_peak_call && !has_output_of_peak_call &&
			 peak_caller_ == 'macs2' && !align_only ) {
			call macs2 { input :
				tas = flatten([[ta_[i]], chosen_ctl_tas[i]]),
				gensz = gensz_,
				chrsz = chrsz_,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				fraglen = fraglen_tmp[i],
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
		File? peak_ = if has_output_of_peak_call then peaks[i]
			else if peak_caller_ == 'spp' then spp.rpeak
			else macs2.npeak

		# signal track
		if ( has_input_of_peak_call && !align_only ) {
			call macs2_signal_track { input :
				tas = flatten([[ta_[i]], chosen_ctl_tas[i]]),
				gensz = gensz_,
				chrsz = chrsz_,
				pval_thresh = pval_thresh,
				fraglen = fraglen_tmp[i],

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}

		# call peaks on 1st pseudo replicated tagalign
		Boolean has_input_of_peak_call_pr1 = defined(spr.ta_pr1[i])
		Boolean has_output_of_peak_call_pr1 = i<length(peaks_pr1) && defined(peaks_pr1[i])
		if ( has_input_of_peak_call_pr1 && !has_output_of_peak_call_pr1 &&
			 peak_caller_=='macs2' && !true_rep_only ) {
			call macs2 as macs2_pr1 { input :
				tas = flatten([[spr.ta_pr1[i]], chosen_ctl_tas[i]]),
				gensz = gensz_,
				chrsz = chrsz_,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				fraglen = fraglen_tmp[i],
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
	
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
		if ( has_input_of_peak_call_pr1 && !has_output_of_peak_call_pr1 &&
			 peak_caller_=='spp' && !true_rep_only ) {
			# call peaks on 1st pseudo replicated tagalign 
			call spp as spp_pr1 { input :
				tas = flatten([[spr.ta_pr1[i]], chosen_ctl_tas[i]]),
				chrsz = chrsz_,
				cap_num_peak = spp_cap_num_peak,
				fraglen = fraglen_tmp[i],
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
	
				cpu = spp_cpu,
				mem_mb = spp_mem_mb,
				disks = spp_disks,
				time_hr = spp_time_hr,
			}
		}
		File? peak_pr1_ = if has_output_of_peak_call_pr1 then peaks_pr1[i]
			else if peak_caller_ == 'spp' then spp_pr1.rpeak
			else macs2_pr1.npeak

		# call peaks on 2nd pseudo replicated tagalign
		Boolean has_input_of_peak_call_pr2 = defined(spr.ta_pr2[i])
		Boolean has_output_of_peak_call_pr2 = i<length(peaks_pr2) && defined(peaks_pr2[i])
		if ( has_input_of_peak_call_pr2 && !has_output_of_peak_call_pr2 &&
			 peak_caller_=='macs2' && !true_rep_only ) {
			call macs2 as macs2_pr2 { input :
				tas = flatten([[spr.ta_pr2[i]], chosen_ctl_tas[i]]),
				gensz = gensz_,
				chrsz = chrsz_,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				fraglen = fraglen_tmp[i],
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
		if ( has_input_of_peak_call_pr2 && !has_output_of_peak_call_pr2 &&
			 peak_caller_=='spp' && !true_rep_only ) {
			# call peaks on 2nd pseudo replicated tagalign 
			call spp as spp_pr2 { input :
				tas = flatten([[spr.ta_pr2[i]], chosen_ctl_tas[i]]),
				chrsz = chrsz_,
				cap_num_peak = spp_cap_num_peak,
				fraglen = fraglen_tmp[i],
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
	
				cpu = spp_cpu,
				mem_mb = spp_mem_mb,
				disks = spp_disks,
				time_hr = spp_time_hr,
			}
		}
		File? peak_pr2_ = if has_output_of_peak_call_pr2 then peaks_pr2[i]
			else if peak_caller_ == 'spp' then spp_pr2.rpeak
			else macs2_pr2.npeak
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

	Boolean has_input_of_peak_caller_pooled = defined(pool_ta.ta_pooled)
	Boolean has_output_of_peak_caller_pooled = defined(peak_pooled)
	if ( has_input_of_peak_caller_pooled && !has_output_of_peak_caller_pooled &&
		peak_caller_=='macs2' && !align_only && num_rep>1 ) {
		# call peaks on pooled replicate
		# always call MACS2 peaks for pooled replicate to get signal tracks
		call macs2 as macs2_pooled { input :
			tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz_,
			chrsz = chrsz_,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( has_input_of_peak_caller_pooled && !has_output_of_peak_caller_pooled &&
		peak_caller_=='spp' && !align_only && num_rep>1 ) {
		# call peaks on pooled replicate
		call spp as spp_pooled { input :
			tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
			chrsz = chrsz_,
			cap_num_peak = spp_cap_num_peak,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			cpu = spp_cpu,
			mem_mb = spp_mem_mb,
			disks = spp_disks,
			time_hr = spp_time_hr,
		}
	}
	File? peak_pooled_ = if has_output_of_peak_caller_pooled then peak_pooled
		else if peak_caller_=='spp' then spp_pooled.rpeak
		else macs2_pooled.npeak	

	# macs2 signal track for pooled rep
	if ( has_input_of_peak_caller_pooled && !align_only && num_rep>1 ) {
		call macs2_signal_track as macs2_signal_track_pooled { input :
			tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz_,
			chrsz = chrsz_,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	Boolean has_input_of_peak_caller_ppr1 = defined(pool_ta_pr1.ta_pooled)
	Boolean has_output_of_peak_caller_ppr1 = defined(peak_ppr1)
	if ( has_input_of_peak_caller_ppr1 && !has_output_of_peak_caller_ppr1 &&
		peak_caller_=='macs2' && !align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 1st pooled pseudo replicates
		call macs2 as macs2_ppr1 { input :
			tas = flatten([select_all([pool_ta_pr1.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz_,
			chrsz = chrsz_,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( has_input_of_peak_caller_ppr1 && !has_output_of_peak_caller_ppr1 &&
		peak_caller_=='spp' && !align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 1st pooled pseudo replicates
		call spp as spp_ppr1 { input :
			tas = flatten([select_all([pool_ta_pr1.ta_pooled]), chosen_ctl_ta_pooled]),
			chrsz = chrsz_,
			cap_num_peak = spp_cap_num_peak,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			cpu = spp_cpu,
			mem_mb = spp_mem_mb,
			disks = spp_disks,
			time_hr = spp_time_hr,
		}
	}
	File? peak_ppr1_ = if has_output_of_peak_caller_ppr1 then peak_ppr1
		else if peak_caller_=='spp' then spp_ppr1.rpeak
		else macs2_ppr1.npeak

	Boolean has_input_of_peak_caller_ppr2 = defined(pool_ta_pr2.ta_pooled)
	Boolean has_output_of_peak_caller_ppr2 = defined(peak_ppr2)
	if ( has_input_of_peak_caller_ppr2 && !has_output_of_peak_caller_ppr2 &&
		peak_caller_=='macs2' && !align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 2nd pooled pseudo replicates
		call macs2 as macs2_ppr2 { input :
			tas = flatten([select_all([pool_ta_pr2.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz_,
			chrsz = chrsz_,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( has_input_of_peak_caller_ppr2 && !has_output_of_peak_caller_ppr2 &&
		peak_caller_=='spp' && !align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 2nd pooled pseudo replicates
		call spp as spp_ppr2 { input :
			tas = flatten([select_all([pool_ta_pr2.ta_pooled]), chosen_ctl_ta_pooled]),
			chrsz = chrsz_,
			cap_num_peak = spp_cap_num_peak,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			cpu = spp_cpu,
			mem_mb = spp_mem_mb,
			disks = spp_disks,
			time_hr = spp_time_hr,
		}
	}
	File? peak_ppr2_ = if has_output_of_peak_caller_ppr2 then peak_ppr2
		else if peak_caller_=='spp' then spp_ppr2.rpeak
		else macs2_ppr2.npeak

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
				prefix = 'rep'+(pair.left+1)+"_rep"+(pair.right+1),
				peak1 = peak_[pair.left],
				peak2 = peak_[pair.right],
				peak_pooled = peak_pooled_,
				fraglen = fraglen_mean.rounded_mean,
				peak_type = peak_type,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
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
				prefix = 'rep'+(pair.left+1)+"_rep"+(pair.right+1),
				peak1 = peak_[pair.left],
				peak2 = peak_[pair.right],
				peak_pooled = peak_pooled_,
				fraglen = fraglen_mean.rounded_mean,
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = pool_ta.ta_pooled,
			}
		}
	}

	# overlap on pseudo-replicates (pr1, pr2) for each true replicate
	scatter( i in range(num_rep) ) {
		if ( !align_only && !true_rep_only ) {
			call overlap as overlap_pr { input :
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peak_pr1_[i],
				peak2 = peak_pr2_[i],
				peak_pooled = peak_[i],
				fraglen = fraglen_[i],
				peak_type = peak_type,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = ta_[i],
			}
		}
	}

	scatter( i in range(num_rep) ) {
		if ( !align_only && !true_rep_only && enable_idr ) {
			# IDR on pseduo replicates
			call idr as idr_pr { input :
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peak_pr1_[i],
				peak2 = peak_pr2_[i],
				peak_pooled = peak_[i],
				fraglen = fraglen_[i],
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = ta_[i],
			}
		}
	}

	if ( !align_only && !true_rep_only && num_rep>1 ) {
		# Naive overlap on pooled pseudo replicates
		call overlap as overlap_ppr { input :
			prefix = "ppr",
			peak1 = peak_ppr1_,
			peak2 = peak_ppr2_,
			peak_pooled = peak_pooled_,
			peak_type = peak_type,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist_,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = pool_ta.ta_pooled,
		}
	}

	if ( !align_only && !true_rep_only && num_rep>1 ) {
		# IDR on pooled pseduo replicates
		call idr as idr_ppr { input :
			prefix = "ppr",
			peak1 = peak_ppr1_,
			peak2 = peak_ppr2_,
			peak_pooled = peak_pooled_,
			idr_thresh = idr_thresh,
			peak_type = peak_type,
			fraglen = fraglen_mean.rounded_mean,
			rank = idr_rank,
			blacklist = blacklist_,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = pool_ta.ta_pooled,
		}
	}

	# reproducibility QC for overlap/IDR peaks
	if ( !align_only && !true_rep_only ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap.bfilt_overlap_peak,
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
			peak_type = peak_type,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	if ( !align_only && !true_rep_only && enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr.bfilt_idr_peak,
			peaks_pr = idr_pr.bfilt_idr_peak,
			peak_ppr = idr_ppr.bfilt_idr_peak,
			peak_type = peak_type,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	# Generate final QC report and JSON
	call qc_report { input :
		pipeline_ver = pipeline_ver,
		title = title,
		description = description,
		genome = basename(select_first([genome_tsv, ref_fa_, chrsz_, 'None'])),
		paired_ends = paired_end_,
		ctl_paired_ends = ctl_paired_end_,
		pipeline_type = pipeline_type,
		peak_caller = peak_caller_,
		macs2_cap_num_peak = macs2_cap_num_peak,
		spp_cap_num_peak = spp_cap_num_peak,		
		idr_thresh = idr_thresh,

		flagstat_qcs = bwa.flagstat_qc,
		nodup_flagstat_qcs = filter.flagstat_qc,
		dup_qcs = filter.dup_qc,
		pbc_qcs = filter.pbc_qc,
		xcor_plots = xcor.plot_png,
		xcor_scores = xcor.score,

		ctl_flagstat_qcs = bwa_ctl.flagstat_qc,
		ctl_nodup_flagstat_qcs = filter_ctl.flagstat_qc,
		ctl_dup_qcs = filter_ctl.dup_qc,
		ctl_pbc_qcs = filter_ctl.pbc_qc,

		jsd_plot = fingerprint.plot,
		jsd_qcs = fingerprint.jsd_qcs,

		frip_macs2_qcs = macs2.frip_qc,
		frip_macs2_qcs_pr1 = macs2_pr1.frip_qc,
		frip_macs2_qcs_pr2 = macs2_pr2.frip_qc,
		frip_macs2_qc_pooled = macs2_pooled.frip_qc,
		frip_macs2_qc_ppr1 = macs2_ppr1.frip_qc,
		frip_macs2_qc_ppr2 = macs2_ppr2.frip_qc,

		frip_spp_qcs = spp.frip_qc,
		frip_spp_qcs_pr1 = spp_pr1.frip_qc,
		frip_spp_qcs_pr2 = spp_pr2.frip_qc,
		frip_spp_qc_pooled = spp_pooled.frip_qc,
		frip_spp_qc_ppr1 = spp_ppr1.frip_qc,
		frip_spp_qc_ppr2 = spp_ppr2.frip_qc,

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
	}

	output {
		File report = qc_report.report
		File qc_json = qc_report.qc_json
		Boolean qc_json_ref_match = qc_report.qc_json_ref_match
	}
}

task merge_fastq { # merge trimmed fastqs
	Array[File] fastqs_R1 		# [merge_id]
	Array[File] fastqs_R2
	Boolean paired_end

	File? null_f
	Array[Array[File]] tmp_fastqs = if paired_end then transpose([fastqs_R1, fastqs_R2])
				else transpose([fastqs_R1])
	command {
		python $(which encode_merge_fastq.py) \
			${write_tsv(tmp_fastqs)} \
			${if paired_end then "--paired-end" else ""} \
			${"--nth " + 1}
	}
	output {
		File merged_fastq_R1 = glob("R1/*.fastq.gz")[0]
		File? merged_fastq_R2 = if paired_end then glob("R2/*.fastq.gz")[0] else null_f
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 2
		disks : "local-disk 100 HDD"
	}
}

task trim_fastq { # trim fastq (for PE R1 only)
	File fastq
	Int trim_bp

	command {
		python $(which encode_trim_fastq.py) \
			${fastq} \
			--trim-bp ${trim_bp}
	}
	output {
		File trimmed_fastq = glob("*.fastq.gz")[0]
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task bwa {
	File bwa_idx_tar	# reference bwa index tar
	File? fastq_R1 		# [read_end_id]
	File? fastq_R2
	Boolean paired_end
	Boolean use_bwa_mem_for_pe

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bwa.py) \
			${bwa_idx_tar} \
			${fastq_R1} ${fastq_R2} \
			${if paired_end then "--paired-end" else ""} \
			${if use_bwa_mem_for_pe then "--use-bwa-mem-for-pe" else ""} \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
		preemptible: 0
	}
}

task filter {
	File bam
	Boolean paired_end
	String dup_marker 				# picard.jar MarkDuplicates (picard) or 
									# sambamba markdup (sambamba)
	Int mapq_thresh				# threshold for low MAPQ reads removal
	Boolean no_dup_removal 		# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
									# filtered bam with dupes	
	String mito_chr_name
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + 0} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if no_dup_removal then "--no-dup-removal" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--nth " + cpu}
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = glob("*.dup.qc")[0]
		File pbc_qc = glob("*.pbc.qc")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task bam2ta {
	File bam
	Boolean paired_end
	String regex_grep_v_ta   	# Perl-style regular expression pattern 
                        		# to remove matching reads from TAGALIGN
	String mito_chr_name 		# mito chromosome name
	Int subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bam2ta.py) \
			${bam} \
			--disable-tn5-shift \
			${if paired_end then "--paired-end" else ""} \
			${if regex_grep_v_ta!="" then "--regex-grep-v-ta '"+regex_grep_v_ta+"'" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task spr { # make two self pseudo replicates
	File ta
	Boolean paired_end

	Int mem_mb

	command {
		python $(which encode_spr.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""}
	}
	output {
		File ta_pr1 = glob("*.pr1.tagAlign.gz")[0]
		File ta_pr2 = glob("*.pr2.tagAlign.gz")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task pool_ta {
	Array[File?] tas

	command {
		python $(which encode_pool_ta.py) \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
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
		python $(which encode_xcor.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
			${"--chip-seq-type " + chip_seq_type} \
			${"--exclusion-range-min " + exclusion_range_min} \
			${"--exclusion-range-max " + exclusion_range_max} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File plot_pdf = glob("*.cc.plot.pdf")[0]
		File plot_png = glob("*.cc.plot.png")[0]
		File score = glob("*.cc.qc")[0]
		Int fraglen = read_int(glob("*.cc.fraglen.txt")[0])
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task fingerprint {
	Array[File?] nodup_bams
	File ctl_bam	 		# one control bam is required
	File blacklist

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_fingerprint.py) \
			${sep=' ' nodup_bams} \
			--ctl-bam ${ctl_bam} \
			${"--blacklist "+ blacklist} \
			${"--nth " + cpu}
	}
	output {
		File plot = glob("*.png")[0]
		Array[File] jsd_qcs = glob("*.jsd.qc")
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
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
		python $(which encode_choose_ctl.py) \
			--tas ${sep=' ' tas} \
			--ctl-tas ${sep=' ' ctl_tas} \
			${"--ta-pooled " + ta_pooled} \
			${"--ctl-ta-pooled " + ctl_ta_pooled} \
			${if always_use_pooled_ctl then "--always-use-pooled-ctl" else ""} \
			${"--ctl-depth-ratio " + ctl_depth_ratio}
	}
	output {
		Array[File] chosen_ctl_tas = glob("ctl_for_rep*.tagAlign.gz")
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task count_signal_track {
	File ta 			# tag-align
	File chrsz			# 2-col chromosome sizes file

	command {
		python $(which encode_count_signal_track.py) \
			${ta} \
			${"--chrsz " + chrsz}
	}
	output {
		File pos_bw = glob("*.positive.bigwig")[0]
		File neg_bw = glob("*.negative.bigwig")[0]
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 4
		disks : "local-disk 50 HDD"
	}
}

task macs2 {
	Array[File?] tas		# [ta, control_ta]. control_ta is optional
	Int fraglen 		# fragment length from xcor
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int cap_num_peak	# cap number of raw peaks called from MACS2
	Float pval_thresh 	# p.value threshold
	File blacklist 		# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak

	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_macs2_chip.py) \
			${sep=' ' tas} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--cap-num-peak " + cap_num_peak} \
			${"--pval-thresh "+ pval_thresh} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--blacklist "+ blacklist}
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t].narrowPeak.gz")[0]
		File bfilt_npeak = glob("*.bfilt.narrowPeak.gz")[0]
		File bfilt_npeak_bb = glob("*.bfilt.narrowPeak.bb")[0]
		File bfilt_npeak_hammock = glob("*.bfilt.narrowPeak.hammock.gz*")[0]
		File bfilt_npeak_hammock_tbi = glob("*.bfilt.narrowPeak.hammock.gz*")[1]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
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
		python $(which encode_macs2_signal_track_chip.py) \
			${sep=' ' tas} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--pval-thresh "+ pval_thresh}
	}
	output {
		File pval_bw = glob("*.pval.signal.bigwig")[0]
		File fc_bw = glob("*.fc.signal.bigwig")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task spp {
	Array[File?] tas		# [ta, control_ta]. control_ta is always required
	Int fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	Int cap_num_peak 	# cap number of raw peaks called from MACS2
	File blacklist 		# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak

	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_spp.py) \
			${sep=' ' tas} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--cap-num-peak " + cap_num_peak} \
			${"--nth " + cpu} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--blacklist "+ blacklist}
	}
	output {
		File rpeak = glob("*[!.][!b][!f][!i][!l][!t].regionPeak.gz")[0]
		File bfilt_rpeak = glob("*.bfilt.regionPeak.gz")[0]
		File bfilt_rpeak_bb = glob("*.bfilt.regionPeak.bb")[0]
		File bfilt_rpeak_hammock = glob("*.bfilt.regionPeak.hammock.gz*")[0]
		File bfilt_rpeak_hammock_tbi = glob("*.bfilt.regionPeak.hammock.gz*")[1]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
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
	File blacklist 		# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
		${if defined(ta) then "" else "touch null.frip.qc"}
		touch null 
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_idr_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		File bfilt_idr_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")[0]
		File bfilt_idr_peak_hammock_tbi = glob("*.bfilt."+peak_type+".hammock.gz*")[1]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.idr*.log")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task overlap {
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	# parameters to compute FRiP
	File? ta		# to calculate FRiP
	Int fraglen 	# fragment length from xcor (for FRIP)
	File chrsz		# 2-col chromosome sizes file
	String peak_type

	command {
		${if defined(ta) then "" else "touch null.frip.qc"}
		touch null 
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			--nonamecheck \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_overlap_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		File bfilt_overlap_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")[0]
		File bfilt_overlap_peak_hammock_tbi = glob("*.bfilt."+peak_type+".hammock.gz*")[1]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
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
	Boolean	keep_irregular_chr_in_bfilt_peak

	command {
		python $(which encode_reproducibility_qc.py) \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix} \
			${"--peak-type " + peak_type} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--chrsz " + chrsz}
	}
	output {
		File optimal_peak = glob("optimal_peak.*.gz")[0]
		File conservative_peak = glob("conservative_peak.*.gz")[0]
		File optimal_peak_bb = glob("optimal_peak.*.bb")[0]
		File conservative_peak_bb = glob("conservative_peak.*.bb")[0]
		File optimal_peak_hammock = glob("optimal_peak.*.hammock.gz*")[0]
		File optimal_peak_hammock_tbi = glob("optimal_peak.*.hammock.gz*")[1]
		File conservative_peak_hammock = glob("conservative_peak.*.hammock.gz*")[0]
		File conservative_peak_hammock_tbi = glob("conservative_peak.*.hammock.gz*")[1]
		File reproducibility_qc = glob("*reproducibility.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
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
	String peak_caller
	Int? macs2_cap_num_peak
	Int? spp_cap_num_peak
	Float idr_thresh
	# QCs
	Array[File?] flagstat_qcs
	Array[File?] nodup_flagstat_qcs
	Array[File?] dup_qcs
	Array[File?] pbc_qcs
	Array[File?] ctl_flagstat_qcs
	Array[File?] ctl_nodup_flagstat_qcs
	Array[File?] ctl_dup_qcs
	Array[File?] ctl_pbc_qcs
	Array[File?] xcor_plots
	Array[File?] xcor_scores
	File? jsd_plot 
	Array[File?] jsd_qcs
	Array[File]? idr_plots
	Array[File?] idr_plots_pr
	File? idr_plot_ppr
	Array[File?] frip_macs2_qcs
	Array[File?] frip_macs2_qcs_pr1
	Array[File?] frip_macs2_qcs_pr2
	File? frip_macs2_qc_pooled
	File? frip_macs2_qc_ppr1 
	File? frip_macs2_qc_ppr2 
	Array[File?] frip_spp_qcs
	Array[File?] frip_spp_qcs_pr1
	Array[File?] frip_spp_qcs_pr2
	File? frip_spp_qc_pooled
	File? frip_spp_qc_ppr1 
	File? frip_spp_qc_ppr2 
	Array[File]? frip_idr_qcs
	Array[File?] frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File?] frip_overlap_qcs
	Array[File?] frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	File? qc_json_ref

	command {
		python $(which encode_qc_report.py) \
			${"--pipeline-ver " + pipeline_ver} \
			${"--title '" + sub(title,"'","_") + "'"} \
			${"--desc '" + sub(description,"'","_") + "'"} \
			${"--genome " + genome} \
			${"--multimapping " + 0} \
			--paired-ends ${sep=" " paired_ends} \
			--pipeline-type ${pipeline_type} \
			--peak-caller ${peak_caller} \
			${"--macs2-cap-num-peak " + macs2_cap_num_peak} \
			${"--spp-cap-num-peak " + spp_cap_num_peak} \
			--idr-thresh ${idr_thresh} \
			--flagstat-qcs ${sep="_:_" flagstat_qcs} \
			--nodup-flagstat-qcs ${sep="_:_" nodup_flagstat_qcs} \
			--dup-qcs ${sep="_:_" dup_qcs} \
			--pbc-qcs ${sep="_:_" pbc_qcs} \
			--xcor-plots ${sep="_:_" xcor_plots} \
			--xcor-scores ${sep="_:_" xcor_scores} \
			--idr-plots ${sep="_:_" idr_plots} \
			--idr-plots-pr ${sep="_:_" idr_plots_pr} \
			--ctl-flagstat-qcs ${sep='_:_' ctl_flagstat_qcs} \
			--ctl-nodup-flagstat-qcs ${sep='_:_' ctl_nodup_flagstat_qcs} \
			--ctl-dup-qcs ${sep='_:_' ctl_dup_qcs} \
			--ctl-pbc-qcs ${sep='_:_' ctl_pbc_qcs} \
			${"--jsd-plot " + jsd_plot} \
			--jsd-qcs ${sep='_:_' jsd_qcs} \
			--frip-spp-qcs ${sep='_:_' frip_spp_qcs} \
			--frip-spp-qcs-pr1 ${sep='_:_' frip_spp_qcs_pr1} \
			--frip-spp-qcs-pr2 ${sep='_:_' frip_spp_qcs_pr2} \
			${"--frip-spp-qc-pooled " + frip_spp_qc_pooled} \
			${"--frip-spp-qc-ppr1 " + frip_spp_qc_ppr1} \
			${"--frip-spp-qc-ppr2 " + frip_spp_qc_ppr2} \
			${"--idr-plot-ppr " + idr_plot_ppr} \
			--frip-macs2-qcs ${sep="_:_" frip_macs2_qcs} \
			--frip-macs2-qcs-pr1 ${sep="_:_" frip_macs2_qcs_pr1} \
			--frip-macs2-qcs-pr2 ${sep="_:_" frip_macs2_qcs_pr2} \
			${"--frip-macs2-qc-pooled " + frip_macs2_qc_pooled} \
			${"--frip-macs2-qc-ppr1 " + frip_macs2_qc_ppr1} \
			${"--frip-macs2-qc-ppr2 " + frip_macs2_qc_ppr2} \
			--frip-idr-qcs ${sep="_:_" frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep="_:_" frip_idr_qcs_pr} \
			${"--frip-idr-qc-ppr " + frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep="_:_" frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep="_:_" frip_overlap_qcs_pr} \
			${"--frip-overlap-qc-ppr " + frip_overlap_qc_ppr} \
			${"--idr-reproducibility-qc " + idr_reproducibility_qc} \
			${"--overlap-reproducibility-qc " + overlap_reproducibility_qc} \
			--out-qc-html qc.html \
			--out-qc-json qc.json \
			${"--qc-json-ref " + qc_json_ref}
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		Boolean qc_json_ref_match = read_string("qc_json_ref_match.txt")=="True"
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

### workflow system tasks

task read_genome_tsv {
	File genome_tsv

	String? null_s
	command <<<
		# create empty files for all entries
		touch ref_fa bowtie2_idx_tar chrsz gensz blacklist

		python <<CODE
		import os
		with open("${genome_tsv}",'r') as fp:
			for line in fp:
				arr = line.strip('\n').split('\t')
				if arr:
					key, val = arr
					with open(key,'w') as fp2:
						fp2.write(val)
		CODE
	>>>
	output {
		String? ref_fa = if size('ref_fa')==0 then null_s else read_string('ref_fa')
		String? bwa_idx_tar = if size('bwa_idx_tar')==0 then null_s else read_string('bwa_idx_tar')
		String? chrsz = if size('chrsz')==0 then null_s else read_string('chrsz')
		String? gensz = if size('gensz')==0 then null_s else read_string('gensz')
		String? blacklist = if size('blacklist')==0 then null_s else read_string('blacklist')
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
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
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task compare_md5sum {
	Array[String] labels
	Array[File] files
	Array[File] ref_files

	command <<<
		python <<CODE	
		from collections import OrderedDict
		import os
		import json
		import hashlib

		def md5sum(filename, blocksize=65536):
		    hash = hashlib.md5()
		    with open(filename, 'rb') as f:
		        for block in iter(lambda: f.read(blocksize), b""):
		            hash.update(block)
		    return hash.hexdigest()

		with open('${write_lines(labels)}','r') as fp:
			labels = fp.read().splitlines()
		with open('${write_lines(files)}','r') as fp:
			files = fp.read().splitlines()
		with open('${write_lines(ref_files)}','r') as fp:
			ref_files = fp.read().splitlines()

		result = OrderedDict()
		match = OrderedDict()
		match_overall = True

		result['tasks'] = []
		result['failed_task_labels'] = []
		result['succeeded_task_labels'] = []
		for i, label in enumerate(labels):
			f = files[i]
			ref_f = ref_files[i]
			md5 = md5sum(f)
			ref_md5 = md5sum(ref_f)
			# if text file, read in contents
			if f.endswith('.qc') or f.endswith('.txt') or \
				f.endswith('.log') or f.endswith('.out'):
				with open(f,'r') as fp:
					contents = fp.read()
				with open(ref_f,'r') as fp:
					ref_contents = fp.read()
			else:
				contents = ''
				ref_contents = ''
			matched = md5==ref_md5
			result['tasks'].append(OrderedDict([
				('label', label),
				('match', matched),
				('md5sum', md5),
				('ref_md5sum', ref_md5),
				('basename', os.path.basename(f)),
				('ref_basename', os.path.basename(ref_f)),
				('contents', contents),
				('ref_contents', ref_contents),
				]))
			match[label] = matched
			match_overall &= matched
			if matched:
				result['succeeded_task_labels'].append(label)
			else:
				result['failed_task_labels'].append(label)		
		result['match_overall'] = match_overall

		with open('result.json','w') as fp:
			fp.write(json.dumps(result, indent=4))
		match_tmp = []
		for key in match:
			val = match[key]
			match_tmp.append('{}\t{}'.format(key, val))
		with open('match.tsv','w') as fp:
			fp.writelines('\n'.join(match_tmp))
		with open('match_overall.txt','w') as fp:
			fp.write(str(match_overall))
		CODE
	>>>
	output {
		Map[String,String] match = read_map('match.tsv') # key:label, val:match
		Boolean match_overall = read_boolean('match_overall.txt')
		File json = glob('result.json')[0] # details (json file)
		String json_str = read_string('result.json') # details (string)
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}
