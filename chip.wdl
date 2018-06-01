# ENCODE DCC TF/Histone ChIP-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow chip {
	#### input file definition
	# pipeline can start from any type of inputs and then leave all other types undefined
	# supported types: fastq, bam, nodup_bam (filtered bam), ta (tagAlign), peak
	# define up to 4 replicates
	# [rep_id] is for each replicate

 	### fastqs
 	# define fastqs either with DNANexus style (1-dim array) or with default one (3-dim array)
 	# [merge_id] is for pooing fastqs
 	## DNANexus UI style fastq definition
	Array[File] fastqs_rep1_R1 = []	# [merge_id]
	Array[File] fastqs_rep1_R2 = [] # do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep2_R1 = [] # do not define if you have a single replicate
	Array[File] fastqs_rep2_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep3_R1 = [] # do not define if you have <=2 replicates
	Array[File] fastqs_rep3_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep4_R1 = [] # do not define if you have <=3 replicates
	Array[File] fastqs_rep4_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] ctl_fastqs_rep1_R1 = []	# [merge_id]
	Array[File] ctl_fastqs_rep1_R2 = [] # do not define _R2 array if your sample is not paired end
	Array[File] ctl_fastqs_rep2_R1 = [] # do not define if you have a single replicate
	Array[File] ctl_fastqs_rep2_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] ctl_fastqs_rep3_R1 = [] # do not define if you have <=2 replicates
	Array[File] ctl_fastqs_rep3_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] ctl_fastqs_rep4_R1 = [] # do not define if you have <=3 replicates
	Array[File] ctl_fastqs_rep4_R2 = []	# do not define _R2 array if your sample is not paired end
 	## default style fastq definition
 	# [read_end_id] is for fastq R1 or fastq R2
	Array[Array[Array[File]]] fastqs = [] 		# [rep_id][merge_id][read_end_id]
	Array[Array[Array[File]]] ctl_fastqs = [] 	# [rep_id][merge_id][read_end_id]

	### other input types (bam, nodup_bam, ta)
	Array[File] bams = [] 			# [rep_id]
	Array[File] ctl_bams = [] 		# [rep_id]
	Array[File] nodup_bams = [] 	# [rep_id]
	Array[File] ctl_nodup_bams = [] # [rep_id]
	Array[File] tas = []			# [rep_id]
	Array[File] ctl_tas = []		# [rep_id]

	### other input types (peak)
	Array[Int] fraglen = [] 	# [rep_id]. fragment length if inputs are peaks
	Array[File] peaks = []		# [PAIR(rep_id1,rep_id2)]. example for 3 reps: [rep1_rep2, rep1_rep3, rep2_rep3]
	Array[File] peaks_pr1 = []	# [rep_id]. do not define if true_rep=true
	Array[File] peaks_pr2 = []	# [rep_id]. do not define if true_rep=true
	File? peak_ppr1				# do not define if you have a single replicate or true_rep=true
	File? peak_ppr2				# do not define if you have a single replicate or true_rep=true
	File? peak_pooled			# do not define if you have a single replicate or true_rep=true

	### pipeline type
	String pipeline_type  		# tf or histone chip-eq
	String peak_caller = '' 	# default: spp for tf and macs2 for histone 

	### mandatory genome param
	File genome_tsv 			# reference genome data TSV file including
								# all important genome specific data file paths and parameters
	Boolean paired_end

	### optional but important
	Boolean align_only = false 		# disable all post-align analysis (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses for pseudo replicates
									# overlap and idr will also be disabled
	Boolean no_jsd = false 			# no JSD plot generation (deeptools fingerprint)

	### task-specific variables but defined in workflow level (limit of WDL)
	## optional for MACS2
	Int macs2_cap_num_peak = 500000	# cap number of raw peaks called from MACS2
	Float pval_thresh = 0.01		# p.value threshold
	Int? macs2_mem_mb 				# resource (memory in MB)
	Int? macs2_time_hr				# resource (walltime in hour)
	String? macs2_disks 			# resource disks for cloud platforms
	## optional for IDR
	Float idr_thresh = 0.05			# IDR threshold
	## optional for SPP
	Int spp_cap_num_peak = 300000	# cap number of raw peaks called from SPP
	Int? spp_cpu	 				# resource (cpu)
	Int? spp_mem_mb 				# resource (memory in MB)
	Int? spp_time_hr				# resource (walltime in hour)
	String? spp_disks 				# resource disks for cloud platforms

	### temp vars (do not define these)
	String peak_caller_ = if peak_caller!='' then peak_caller 
						else if pipeline_type=='tf' then 'spp'
						else 'macs2'
	String peak_type = if peak_caller_=='spp' then 'regionPeak'
						else if peak_caller_=='macs2' then 'narrowPeak'
						else 'narrowPeak'
	Boolean enable_idr = pipeline_type=='tf' # enable_idr for TF chipseq only
	String idr_rank = if peak_caller_=='spp' then 'signal.value'
						else if peak_caller_=='macs2' then 'p.value'
						else 'p.value'

	### read genome data and paths
	call read_genome_tsv { input:genome_tsv = genome_tsv }
	File bwa_idx_tar = read_genome_tsv.genome['bwa_idx_tar']
	File blacklist = read_genome_tsv.genome['blacklist']
	File chrsz = read_genome_tsv.genome['chrsz']
	String gensz = read_genome_tsv.genome['gensz']

	### pipeline starts here
	# temporary 2-dim arrays for DNANexus style fastqs and adapters
	Array[Array[File]] fastqs_rep1 = transpose([fastqs_rep1_R1,fastqs_rep1_R2])
	Array[Array[File]] fastqs_rep2 = transpose([fastqs_rep2_R1,fastqs_rep2_R2])
	Array[Array[File]] fastqs_rep3 = transpose([fastqs_rep3_R1,fastqs_rep3_R2])
	Array[Array[File]] fastqs_rep4 = transpose([fastqs_rep4_R1,fastqs_rep4_R2])

	Array[Array[Array[File]]] fastqs_ = if length(fastqs_rep1)<1 then fastqs
		else if length(fastqs_rep2)<1 then [fastqs_rep1]
		else if length(fastqs_rep3)<1 then [fastqs_rep1,fastqs_rep2]
		else if length(fastqs_rep4)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3]
		else [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4]
	scatter(fastq_set in fastqs_) {
		# merge fastqs
		call merge_fastq { input :
			fastqs = fastq_set,
			paired_end = paired_end,
		}
		# align merged fastqs with bwa
		call bwa { input :
			idx_tar = bwa_idx_tar,
			fastqs = merge_fastq.merged_fastqs, #[R1,R2]
			paired_end = paired_end,
		}
	}

	# special treatment for xcor for paired end samples only
	Array[Array[File]] fastqs_xcor = if !paired_end then [] else merge_fastq.merged_fastqs
	scatter(fastq_set in fastqs_xcor) {
		# for paired end dataset, map R1 only as SE for xcor analysis
		call trim_fastq { input :
			fastq = fastq_set[0],
		}
	}
	Array[Array[File]] trimmed_fastqs_R1 = if length(trim_fastq.trimmed_fastq)<1 then []
		else transpose([trim_fastq.trimmed_fastq])
	scatter(fastq_set in trimmed_fastqs_R1) {
		call bwa as bwa_R1 { input :
			idx_tar = bwa_idx_tar,
			fastqs = fastq_set,
			paired_end = false,
		}
		# no bam filtering for xcor
		call bam2ta as bam2ta_no_filt_R1 { input :
			bam = bwa_R1.bam,
			disable_tn5_shift = true,
			paired_end = false,
		}
	}

	Array[File] bams_ = flatten([bwa.bam, bams])
	scatter(bam in bams_) {
		# filter/dedup bam
		call filter { input :
			bam = bam,
			paired_end = paired_end,
		}
		# make unfiltered bam for xcor
		call bam2ta as bam2ta_no_filt { input :
			bam = bam,
			disable_tn5_shift = true,
			paired_end = paired_end,
		}
	}

	Array[File] nodup_bams_ = flatten([filter.nodup_bam, nodup_bams])
	scatter(bam in nodup_bams_) {
		# convert bam to tagalign and subsample it if necessary
		call bam2ta { input :
			bam = bam,
			disable_tn5_shift = true,
			paired_end = paired_end,
		}
	}

	Array[File] tas_xcor = if length(bam2ta_no_filt_R1.ta)>0 then bam2ta_no_filt_R1.ta
		else if length(bam2ta_no_filt.ta)>0 then bam2ta_no_filt.ta
		else flatten([bam2ta.ta, tas])
	Boolean paired_end_xcor = paired_end && length(bam2ta_no_filt_R1.ta)<1
	scatter(ta in tas_xcor) {
		# use trimmed/unfilitered R1 tagAlign for paired end dataset 		
		# if not starting from fastqs, keep using old method
		#  (mapping with both ends for tag-aligns to be used for xcor)
		# subsample tagalign (non-mito) and cross-correlation analysis
		call xcor { input :
			ta = ta,
			paired_end = paired_end_xcor,
		}
	}

	Array[File] tas_ = if align_only then [] else flatten([bam2ta.ta, tas])	
	if ( length(tas_)>1 )  {
		# pool tagaligns from true replicates
		call pool_ta { input :
			tas = tas_,
		}
	}

	if ( !true_rep_only ) {
		scatter( ta in tas_ ) {
			# make two self pseudo replicates per true replicate
			call spr { input :
				ta = ta,
				paired_end = paired_end,
			}
		}
	}
	if ( !true_rep_only && length(tas_)>1 ) {
		# pool tagaligns from pseudo replicates
		call pool_ta as pool_ta_pr1 { input :
			tas = spr.ta_pr1,
		}
		call pool_ta as pool_ta_pr2 { input :
			tas = spr.ta_pr2,
		}
	}

	# align controls
	Array[Array[File]] ctl_fastqs_rep1 = transpose([ctl_fastqs_rep1_R1,ctl_fastqs_rep1_R2])
	Array[Array[File]] ctl_fastqs_rep2 = transpose([ctl_fastqs_rep2_R1,ctl_fastqs_rep2_R2])
	Array[Array[File]] ctl_fastqs_rep3 = transpose([ctl_fastqs_rep3_R1,ctl_fastqs_rep3_R2])
	Array[Array[File]] ctl_fastqs_rep4 = transpose([ctl_fastqs_rep4_R1,ctl_fastqs_rep4_R2])
	Array[Array[Array[File]]] ctl_fastqs_ = if length(ctl_fastqs_rep1)<1 then ctl_fastqs
		else if length(ctl_fastqs_rep2)<1 then [ctl_fastqs_rep1]
		else if length(ctl_fastqs_rep3)<1 then [ctl_fastqs_rep1,ctl_fastqs_rep2]
		else if length(ctl_fastqs_rep4)<1 then [ctl_fastqs_rep1,ctl_fastqs_rep2,ctl_fastqs_rep3]
		else [ctl_fastqs_rep1,ctl_fastqs_rep2,ctl_fastqs_rep3,ctl_fastqs_rep4]
	scatter(fastq_set in ctl_fastqs_) {
		# merge fastqs
		call merge_fastq as merge_fastq_ctl { input :
			fastqs = fastq_set,
			paired_end = paired_end,
		}
		# align merged fastqs with bwa
		call bwa as bwa_ctl { input :
			idx_tar = bwa_idx_tar,
			fastqs = merge_fastq_ctl.merged_fastqs, #[R1,R2]
			paired_end = paired_end,
		}
	}

	Array[File] ctl_bams_ = flatten([bwa_ctl.bam, ctl_bams])
	scatter(bam in ctl_bams_) {
		# filter/dedup bam
		call filter as filter_ctl { input :
			bam = bam,
			paired_end = paired_end,
		}
	}

	Array[File] ctl_nodup_bams_ = flatten([filter_ctl.nodup_bam, ctl_nodup_bams])
	scatter(bam in ctl_nodup_bams_) {
		# convert bam to tagalign and subsample it if necessary
		call bam2ta as bam2ta_ctl { input :
			bam = bam,
			disable_tn5_shift = true,
			paired_end = paired_end,
		}
	}

	Array[String] ctl_tas_ = if align_only then [] else flatten([bam2ta_ctl.ta, ctl_tas])
	if ( length(ctl_tas_)>0 ) {	
		# pool tagaligns from true replicates
		call pool_ta as pool_ta_ctl { input :
			tas = ctl_tas_,
		}
	}

	if ( !no_jsd && length(nodup_bams_)>0 && length(ctl_nodup_bams_)>0 && basename(blacklist)!='null' ) {
		# fingerprint and JS-distance plot
		call fingerprint { input :
			nodup_bams = nodup_bams_,
			ctl_bam = ctl_nodup_bams_[0], # use first control only
			blacklist = blacklist,
		}
	}

	if ( length(tas_)>0 && length(ctl_tas_)>0 ) {
		# choose appropriate control for each exp IP replicate
		# outputs:
		# 	choose_ctl.idx : control replicate index for each exp replicate 
		#					-1 means pooled ctl replicate
		call choose_ctl { input:
			tas = tas_,
			ctl_tas = ctl_tas_,
			ta_pooled = pool_ta.ta_pooled,
			ctl_ta_pooled = pool_ta_ctl.ta_pooled,
		}
	}
	# before peak calling, get fragment length from xcor analysis or given input
	# if fraglen [] is defined in the input JSON, fraglen from xcor will be ignored
	Array[Int] fraglen_ = if align_only then [] 
		else if length(fraglen)>0 then fraglen
		else xcor.fraglen

	# make control ta array [[1,2,3,4]] -> [[1],[2],[3],[4]], will be zipped with exp ta array latter
	Array[Array[File]] chosen_ctl_tas =	if length(tas_)<1 || length(ctl_tas_)<1 then [[],[],[],[],[],[]]
		else transpose(select_all([choose_ctl.chosen_ctl_tas]))

	# we have all tas and ctl_tas (optional for histone chipseq) ready, let's call peaks
	scatter(i in range(length(tas_))) {
		# always call MACS2 peaks for true replicates to get signal tracks
		# call peaks on tagalign
		call macs2 { input :
			tas = flatten([[tas_[i]], chosen_ctl_tas[i]]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			make_signal = true,
			fraglen = fraglen_[i],
			blacklist = blacklist,
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	# SPP cannot call peaks without controls
	if ( peak_caller_=='spp' ) {
		scatter(i in range(length(tas_))) {
			# call peaks on tagalign
			call spp { input :
				tas = flatten([[tas_[i]], chosen_ctl_tas[i]]),
				chrsz = chrsz,
				cap_num_peak = spp_cap_num_peak,
				fraglen = fraglen_[i],
				blacklist = blacklist,
				# resource
				cpu = spp_cpu,
				mem_mb = spp_mem_mb,
				disks = spp_disks,
				time_hr = spp_time_hr,
			}
		}
	}

	if ( peak_caller_=='macs2' ) {
		scatter(i in range(length(tas_))) {
			# call peaks on 1st pseudo replicated tagalign 
			call macs2 as macs2_pr1 { input :
				tas = flatten([[select_first([spr.ta_pr1])[i]], chosen_ctl_tas[i]]),
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				fraglen = fraglen_[i],
				blacklist = blacklist,
				# resource
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
			call macs2 as macs2_pr2 { input :
				tas = flatten([[select_first([spr.ta_pr2])[i]], chosen_ctl_tas[i]]),
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				fraglen = fraglen_[i],
				blacklist = blacklist,
				# resource
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}					
		}
	}

	if ( peak_caller_=='spp' ) {
		scatter(i in range(length(tas_))) {
			# call peaks on 1st pseudo replicated tagalign 
			call spp as spp_pr1 { input :
				tas = flatten([[select_first([spr.ta_pr1])[i]], chosen_ctl_tas[i]]),
				chrsz = chrsz,
				cap_num_peak = spp_cap_num_peak,
				fraglen = fraglen_[i],
				blacklist = blacklist,
				# resource
				cpu = spp_cpu,
				mem_mb = spp_mem_mb,
				disks = spp_disks,
				time_hr = spp_time_hr,
			}
			# call peaks on 2nd pseudo replicated tagalign 
			call spp as spp_pr2 { input :
				tas = flatten([[select_first([spr.ta_pr2])[i]], chosen_ctl_tas[i]]),
				chrsz = chrsz,
				cap_num_peak = spp_cap_num_peak,
				fraglen = fraglen_[i],
				blacklist = blacklist,
				# resource
				cpu = spp_cpu,
				mem_mb = spp_mem_mb,
				disks = spp_disks,
				time_hr = spp_time_hr,
			}
		}
	}

	# rounded mean of fragment length, which will be used for 
	#  1) calling peaks for pooled true/pseudo replicates
	#  2) calculating FRiP
	call rounded_mean as fraglen_mean { input :
		ints = fraglen_,
	}

	# actually not an array
	Array[File] chosen_ctl_ta_pooled = if length(tas_)<2 || length(ctl_tas_)<1 then []
		else if length(ctl_tas_)<2 then [ctl_tas_[0]] # choose first (only) control
		else select_all([pool_ta_ctl.ta_pooled]) # choose pooled control

	if ( length(tas_)>1 ) {
		# call peaks on pooled replicate
		# always call MACS2 peaks for pooled replicate to get signal tracks
		call macs2 as macs2_pooled { input :
			tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			make_signal = true,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( length(tas_)>1 && peak_caller_=='spp' ) {
		# call peaks on pooled replicate
		call spp as spp_pooled { input :
			tas = flatten([select_all([pool_ta.ta_pooled]), chosen_ctl_ta_pooled]),
			chrsz = chrsz,
			cap_num_peak = spp_cap_num_peak,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			# resource
			cpu = spp_cpu,
			mem_mb = spp_mem_mb,
			disks = spp_disks,
			time_hr = spp_time_hr,
		}
	}

	if ( !true_rep_only && length(tas_)>1 && peak_caller_=='macs2' ) {
		# call peaks on 1st pooled pseudo replicates
		call macs2 as macs2_ppr1 { input :
			tas = flatten([select_all([pool_ta_pr1.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
		# call peaks on 2nd pooled pseudo replicates
		call macs2 as macs2_ppr2 { input :
			tas = flatten([select_all([pool_ta_pr2.ta_pooled]), chosen_ctl_ta_pooled]),
			gensz = gensz,
			chrsz = chrsz,
			cap_num_peak = macs2_cap_num_peak,
			pval_thresh = pval_thresh,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( !true_rep_only && length(tas_)>1 && peak_caller_=='spp' ) {
		# call peaks on 1st pooled pseudo replicates
		call spp as spp_ppr1 { input :
			tas = flatten([select_all([pool_ta_pr1.ta_pooled]), chosen_ctl_ta_pooled]),
			chrsz = chrsz,
			cap_num_peak = spp_cap_num_peak,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			# resource
			cpu = spp_cpu,
			mem_mb = spp_mem_mb,
			disks = spp_disks,
			time_hr = spp_time_hr,
		}
		# call peaks on 2nd pooled pseudo replicates
		call spp as spp_ppr2 { input :
			tas = flatten([select_all([pool_ta_pr2.ta_pooled]), chosen_ctl_ta_pooled]),
			chrsz = chrsz,
			cap_num_peak = spp_cap_num_peak,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			# resource
			cpu = spp_cpu,
			mem_mb = spp_mem_mb,
			disks = spp_disks,
			time_hr = spp_time_hr,
		}
	}

	# make peak arrays
	Array[File] peaks_ = if align_only then [] else select_first([spp.rpeak ,flatten([macs2.npeak, peaks])])

	# generate all possible pairs of true replicates (pair: left=prefix, right=[peak1,peak2])
	Array[Pair[String,Array[File]]] peak_pairs =  
		if length(peaks_)<=1 then [] # 1 rep
		else if length(peaks_)<=2 then # 2 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]])]
		else if length(peaks_)<=3 then # 3 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]])]
		else if length(peaks_)<=4 then # 4 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]])]
		else if length(peaks_)<=5 then # 5 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]), ('rep1-rep5',[peaks_[0],peaks_[4]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]), ('rep2-rep5',[peaks_[1],peaks_[4]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]]), ('rep3-rep5',[peaks_[2],peaks_[4]]),
			  ('rep4-rep5',[peaks_[3],peaks_[4]])]
		else # 6 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]), ('rep1-rep5',[peaks_[0],peaks_[4]]), ('rep1-rep6',[peaks_[0],peaks_[5]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]), ('rep2-rep5',[peaks_[1],peaks_[4]]), ('rep2-rep6',[peaks_[1],peaks_[5]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]]), ('rep3-rep5',[peaks_[2],peaks_[4]]), ('rep3-rep6',[peaks_[2],peaks_[5]]),
			  ('rep4-rep5',[peaks_[3],peaks_[4]]), ('rep4-rep6',[peaks_[3],peaks_[5]]),
			  ('rep5-rep6',[peaks_[4],peaks_[5]])]
	scatter( pair in peak_pairs ) {
		# Naive overlap on every pair of true replicates
		call overlap { input :
			prefix = pair.left,
			peak1 = pair.right[0],
			peak2 = pair.right[1],
			peak_pooled = select_first([spp_pooled.rpeak, macs2_pooled.npeak, peak_pooled]),
			fraglen = fraglen_mean.rounded_mean,
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	if ( enable_idr ) {
		scatter( pair in peak_pairs ) {
			# IDR on every pair of true replicates
			call idr { input : 
				prefix = pair.left,
				peak1 = pair.right[0],
				peak2 = pair.right[1],
				peak_pooled = select_first([spp_pooled.rpeak, macs2_pooled.npeak, peak_pooled]),
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				fraglen = fraglen_mean.rounded_mean,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				ta = pool_ta.ta_pooled,
			}
		}
	}

	Array[File] peaks_pr1_ = select_first([spp_pr1.rpeak, macs2_pr1.npeak, peaks_pr1])
	Array[File] peaks_pr2_ = select_first([spp_pr2.rpeak, macs2_pr2.npeak, peaks_pr2])
	scatter( i in range(length(peaks_pr1_)) ) {
		# Naive overlap on pseduo replicates
		call overlap as overlap_pr { input : 
			prefix = "rep"+(i+1)+"-pr",
			peak1 = peaks_pr1_[i],
			peak2 = peaks_pr2_[i],
			peak_pooled = peaks_[i],
			fraglen = fraglen_[i],
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = if length(tas_)>0 then tas_[i] else pool_ta.ta_pooled,
		}
	}
	if ( enable_idr ) {
		scatter( i in range(length(peaks_pr1_)) ) {
			# IDR on pseduo replicates
			call idr as idr_pr { input : 
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peaks_pr1_[i],
				peak2 = peaks_pr2_[i],
				peak_pooled = peaks_[i],
				fraglen = fraglen_[i],
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				ta = if length(tas_)>0 then tas_[i] else pool_ta.ta_pooled,
			}
		}
	}

	if ( length(peaks_pr1_)>1 ) {
		# Naive overlap on pooled pseudo replicates
		call overlap as overlap_ppr { input : 
			prefix = "ppr",
			peak1 = select_first([spp_ppr1.rpeak, macs2_ppr1.npeak, peak_ppr1]), #peak_ppr1_[0],
			peak2 = select_first([spp_ppr2.rpeak, macs2_ppr2.npeak, peak_ppr2]), #peak_ppr2_[0],
			peak_pooled = select_first([spp_pooled.rpeak, macs2_pooled.npeak, peak_pooled]),
			peak_type = peak_type,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	if ( enable_idr && length(peaks_pr1_)>1 ) {
		# IDR on pooled pseduo replicates
		call idr as idr_ppr { input : 
			prefix = "ppr",
			peak1 = select_first([spp_ppr1.rpeak, macs2_ppr1.npeak, peak_ppr1]), #peak_ppr1_[0],
			peak2 = select_first([spp_ppr2.rpeak, macs2_ppr2.npeak, peak_ppr2]), #peak_ppr2_[0],
			peak_pooled = select_first([spp_pooled.rpeak, macs2_pooled.npeak, peak_pooled]),
			idr_thresh = idr_thresh,
			peak_type = peak_type,
			rank = idr_rank,
			fraglen = fraglen_mean.rounded_mean,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}

	if ( !align_only && !true_rep_only ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap.bfilt_overlap_peak,
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
		}
	}
	if ( !align_only && !true_rep_only && enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr.bfilt_idr_peak,
			peaks_pr = idr_pr.bfilt_idr_peak,
			peak_ppr = idr_ppr.bfilt_idr_peak,
		}
	}
	# Generate final QC report and JSON
	call qc_report { input :
		paired_end = paired_end,
		pipeline_type = pipeline_type,
		peak_caller = peak_caller_,
		idr_thresh = idr_thresh,
		flagstat_qcs = bwa.flagstat_qc,
		nodup_flagstat_qcs = filter.flagstat_qc,
		dup_qcs = filter.dup_qc,
		pbc_qcs = filter.pbc_qc,
		ctl_flagstat_qcs = bwa_ctl.flagstat_qc,
		ctl_nodup_flagstat_qcs = filter_ctl.flagstat_qc,
		ctl_dup_qcs = filter_ctl.dup_qc,
		ctl_pbc_qcs = filter_ctl.pbc_qc,
		xcor_plots = xcor.plot_png,
		xcor_scores = xcor.score,

		jsd_plot = fingerprint.plot,
		jsd_qcs = select_first([fingerprint.jsd_qcs,[]]),

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
}

task merge_fastq { # merge trimmed fastqs
	# parameters from workflow
	Array[Array[File]] fastqs 		# [merge_id][read_end_id]
	Boolean paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_merge_fastq.py) \
			${write_tsv(fastqs)} \
			${if paired_end then "--paired-end" else ""} \
			${"--nth " + select_first([cpu,2])}
	}
	output {
		# WDL glob() globs in an alphabetical order
		# so R1 and R2 can be switched, which results in an
		# unexpected behavior of a workflow
		# so we prepend merge_fastqs_'end'_ (R1 or R2)
		# to the basename of original filename
		# this prefix will be later stripped in bwa task
		Array[File] merged_fastqs = glob("merge_fastqs_R?_*.fastq.gz")
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task trim_fastq { # trim fastq (for PE R1 only)
	# parameters from workflow
	File fastq
	Int? trim_bp

	command {
		python $(which encode_trim_fastq.py) \
			${fastq} \
			--trim-bp ${select_first([trim_bp,50])}
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
	# parameters from workflow
	File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [read_end_id]
	Boolean paired_end

	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_bwa.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${if paired_end then "--paired-end" else ""} \
			${"--nth " + select_first([cpu,4])}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : select_first([cpu,4])
		memory : "${select_first([mem_mb,'20000'])} MB"
		time : select_first([time_hr,48])
		disks : select_first([disks,"local-disk 100 HDD"])
		preemptible: 0
	}
}

task filter {
	# parameters from workflow
	File bam
	Boolean paired_end
	Int? multimapping
	# optional
	String? dup_marker 				# picard.jar MarkDuplicates (picard) or 
									# sambamba markdup (sambamba)
	Int? mapq_thresh				# threshold for low MAPQ reads removal
	Boolean? no_dup_removal 		# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
	# resource						# filtered bam with dupes	
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--dup-marker " + select_first([dup_marker,'picard'])} \
			${"--mapq-thresh " + select_first([mapq_thresh,30])} \
			${if select_first([no_dup_removal,false]) then "--no-dup-removal" else ""} \
			${"--nth " + cpu}
		# ugly part to deal with optional outputs with Google JES backend
		${if select_first([no_dup_removal,false]) then "touch null.dup.qc null.pbc.qc; " else ""}
		touch null
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.dup.qc")[0]
		File pbc_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.pbc.qc")[0]
	}
	runtime {
		#@docker : "quay.io/encode-dcc/atac-seq-pipeline:v1"
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'20000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task bam2ta {
	# parameters from workflow
	File bam
	Boolean paired_end
	Boolean disable_tn5_shift 	# no tn5 shifting (it's for dnase-seq)
	# optional
	String? regex_grep_v_ta   	# Perl-style regular expression pattern 
                        		# to remove matching reads from TAGALIGN
	Int? subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_bam2ta.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${if disable_tn5_shift then "--disable-tn5-shift" else ""} \
			${"--regex-grep-v-ta " +"'"+select_first([regex_grep_v_ta,'chrM'])+"'"} \
			${"--subsample " + select_first([subsample,0])} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}
	runtime {
		#@docker : "quay.io/encode-dcc/atac-seq-pipeline:v1"
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task spr { # make two self pseudo replicates
	# parameters from workflow
	File ta
	Boolean paired_end

	# resource
	Int? mem_mb

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
		memory : "${select_first([mem_mb,'12000'])} MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task pool_ta {
	# parameters from workflow
	Array[File] tas

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
	# parameters from workflow
	File ta
	Boolean paired_end
	# optional
	Int? subsample  # number of reads to subsample TAGALIGN
					# this will be used for xcor only
					# will not affect any downstream analysis
	# resource
	Int? cpu
	Int? mem_mb	
	Int? time_hr
	String? disks

	command {
		python $(which encode_xcor.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""} \
			${"--subsample " + select_first([subsample,15000000])} \
			${"--nth " + select_first([cpu,2])}
	}
	output {
		File plot_pdf = glob("*.cc.plot.pdf")[0]
		File plot_png = glob("*.cc.plot.png")[0]
		File score = glob("*.cc.qc")[0]
		Int fraglen = read_int(glob("*.cc.fraglen.txt")[0])
	}
	runtime {
		#@docker : "quay.io/encode-dcc/atac-seq-pipeline:v1"
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task fingerprint {
	# parameters from workflow
	Array[File] nodup_bams
	File ctl_bam	 		# one control bam is required
	File blacklist

	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_fingerprint.py) \
			${sep=' ' nodup_bams} \
			--ctl-bam ${ctl_bam} \
			${"--blacklist "+ blacklist} \
			${"--nth " + select_first([cpu,2])}
	}
	output {
		File plot = glob("*.png")[0]
		Array[File] jsd_qcs = glob("*.jsd.qc")
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task choose_ctl {
	# parameters from workflow
	Array[File] tas
	Array[File] ctl_tas
	File? ta_pooled
	File? ctl_ta_pooled
	# optional
	Boolean? always_use_pooled_ctl # always use pooled control for all exp rep.
	Float? ctl_depth_ratio 		# if ratio between controls is higher than this
								# then always use pooled control for all exp rep.
	command {
		python $(which encode_choose_ctl.py) \
			--tas ${sep=' ' tas} \
			--ctl-tas ${sep=' ' ctl_tas} \
			${"--ta-pooled " + ta_pooled} \
			${"--ctl-ta-pooled " + ctl_ta_pooled} \
			${if select_first([always_use_pooled_ctl,false]) then "--always-use-pooled-ctl" else ""} \
			${"--ctl-depth-ratio " + select_first([ctl_depth_ratio,'1.2'])}
	}
	output {
		#Array[Int] idx = read_lines("idx.txt")
		Array[File] chosen_ctl_tas = glob("ctl_for_rep*.tagAlign.gz")
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task macs2 {
	# parameters from workflow
	Array[File] tas		# [ta, control_ta]. control_ta is optional
	Int fraglen 		# fragment length from xcor
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int? cap_num_peak	# cap number of raw peaks called from MACS2
	Float? pval_thresh 	# p.value threshold
	Boolean? make_signal
	File blacklist 		# blacklist BED to filter raw peaks
	# resource
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_macs2_chip.py) \
			${sep=' ' tas} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--cap-num-peak " + select_first([cap_num_peak,500000])} \
			${"--pval-thresh "+ select_first([pval_thresh,'0.01'])} \
			${if select_first([make_signal,false]) then "--make-signal" else ""} \
			${"--blacklist "+ blacklist}

		${if select_first([make_signal,false]) then "" 
			else "touch null.pval.signal.bigwig null.fc.signal.bigwig"}
		touch null # ugly part to deal with optional outputs
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t].narrowPeak.gz")[0]
		File bfilt_npeak = glob("*.bfilt.narrowPeak.gz")[0]
		File bfilt_npeak_bb = glob("*.bfilt.narrowPeak.bb")[0]
		File sig_pval = if select_first([make_signal,false]) then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if select_first([make_signal,false]) then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task spp {
	# parameters from workflow
	Array[File] tas		# [ta, control_ta]. control_ta is always required
	Int fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	Int? cap_num_peak 	# cap number of raw peaks called from MACS2
	File blacklist 		# blacklist BED to filter raw peaks
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_spp.py) \
			${sep=' ' tas} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--cap-num-peak " + select_first([cap_num_peak,300000])} \
			${"--nth " + select_first([cpu,2])} \
			${"--blacklist "+ blacklist}
	}
	output {
		File rpeak = glob("*[!.][!b][!f][!i][!l][!t].regionPeak.gz")[0]
		File bfilt_rpeak = glob("*.bfilt.regionPeak.gz")[0]
		File bfilt_rpeak_peak_bb = glob("*.bfilt.regionPeak.bb")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,72])
		disks : select_first([disks,"local-disk 100 HDD"])
		preemptible: 0
	}
}

task idr {
	# parameters from workflow
	String? prefix 		# prefix for IDR output file
	File peak1 			
	File peak2
	File peak_pooled
	Float? idr_thresh
	File blacklist 		# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + select_first([idr_thresh,'0.05'])} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(ta) then "" else "touch null.frip.qc"}			
		touch null 
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_idr_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
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
	# parameters from workflow
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File? ta		# to calculate FRiP
	Int fraglen 	# fragment length from xcor (for FRIP)
	File chrsz		# 2-col chromosome sizes file
	String peak_type

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(ta) then "" else "touch null.frip.qc"}			
		touch null 
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_overlap_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
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
	# parameters from workflow
	String prefix
	Array[File] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File] peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.

	command {
		python $(which encode_reproducibility_qc.py) \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix}
	}
	output {
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
 	String? name # name of sample
	String? desc # description for sample
	#String? encode_accession_id	# ENCODE accession ID of sample
	# workflow params
	Boolean paired_end
	String pipeline_type
	String peak_caller
	Float idr_thresh
	# QCs
	Array[File]? flagstat_qcs
	Array[File]? nodup_flagstat_qcs
	Array[File]? dup_qcs
	Array[File]? pbc_qcs
	Array[File]? ctl_flagstat_qcs
	Array[File]? ctl_nodup_flagstat_qcs
	Array[File]? ctl_dup_qcs
	Array[File]? ctl_pbc_qcs
	Array[File]? xcor_plots
	Array[File]? xcor_scores
	File? jsd_plot 
	Array[File]? jsd_qcs
	Array[File]? idr_plots
	Array[File]? idr_plots_pr
	File? idr_plot_ppr
	Array[File]? frip_macs2_qcs
	Array[File]? frip_macs2_qcs_pr1
	Array[File]? frip_macs2_qcs_pr2
	File? frip_macs2_qc_pooled
	File? frip_macs2_qc_ppr1 
	File? frip_macs2_qc_ppr2 
	Array[File]? frip_spp_qcs
	Array[File]? frip_spp_qcs_pr1
	Array[File]? frip_spp_qcs_pr2
	File? frip_spp_qc_pooled
	File? frip_spp_qc_ppr1 
	File? frip_spp_qc_ppr2 
	Array[File]? frip_idr_qcs
	Array[File]? frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File]? frip_overlap_qcs
	Array[File]? frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	File? qc_json_ref

	command {
		python $(which encode_qc_report.py) \
			${"--name '" + name + "'"} \
			${"--desc '" + desc + "'"} \
			${if paired_end then "--paired-end" else ""} \
			--pipeline-type ${pipeline_type} \
			--peak-caller ${peak_caller} \
			--idr-thresh ${idr_thresh} \
			--flagstat-qcs ${sep=' ' flagstat_qcs} \
			--nodup-flagstat-qcs ${sep=' ' nodup_flagstat_qcs} \
			--dup-qcs ${sep=' ' dup_qcs} \
			--pbc-qcs ${sep=' ' pbc_qcs} \
			--ctl-flagstat-qcs ${sep=' ' ctl_flagstat_qcs} \
			--ctl-nodup-flagstat-qcs ${sep=' ' ctl_nodup_flagstat_qcs} \
			--ctl-dup-qcs ${sep=' ' ctl_dup_qcs} \
			--ctl-pbc-qcs ${sep=' ' ctl_pbc_qcs} \
			--xcor-plots ${sep=' ' xcor_plots} \
			--xcor-scores ${sep=' ' xcor_scores} \
			${"--jsd-plot " + jsd_plot} \
			--jsd-qcs ${sep=' ' jsd_qcs} \
			--idr-plots ${sep=' ' idr_plots} \
			--idr-plots-pr ${sep=' ' idr_plots_pr} \
			${"--idr-plot-ppr " + idr_plot_ppr} \
			--frip-macs2-qcs ${sep=' ' frip_macs2_qcs} \
			--frip-macs2-qcs-pr1 ${sep=' ' frip_macs2_qcs_pr1} \
			--frip-macs2-qcs-pr2 ${sep=' ' frip_macs2_qcs_pr2} \
			${"--frip-macs2-qc-pooled " + frip_macs2_qc_pooled} \
			${"--frip-macs2-qc-ppr1 " + frip_macs2_qc_ppr1} \
			${"--frip-macs2-qc-ppr2 " + frip_macs2_qc_ppr2} \
			--frip-spp-qcs ${sep=' ' frip_spp_qcs} \
			--frip-spp-qcs-pr1 ${sep=' ' frip_spp_qcs_pr1} \
			--frip-spp-qcs-pr2 ${sep=' ' frip_spp_qcs_pr2} \
			${"--frip-spp-qc-pooled " + frip_spp_qc_pooled} \
			${"--frip-spp-qc-ppr1 " + frip_spp_qc_ppr1} \
			${"--frip-spp-qc-ppr2 " + frip_spp_qc_ppr2} \
			--frip-idr-qcs ${sep=' ' frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep=' ' frip_idr_qcs_pr} \
			${"--frip-idr-qc-ppr " + frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep=' ' frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep=' ' frip_overlap_qcs_pr} \
			${"--frip-overlap-qc-ppr " + frip_overlap_qc_ppr} \
			${"--idr-reproducibility-qc " + idr_reproducibility_qc} \
			${"--overlap-reproducibility-qc " + overlap_reproducibility_qc} \
			--out-qc-html qc.html \
			--out-qc-json qc.json

		diff qc.json ${if defined(qc_json_ref) then qc_json_ref else "/dev/null"} | wc -l > qc_json_match.txt
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		Boolean qc_json_match = read_int("qc_json_match.txt")==0
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
	command {
		echo "Reading genome_tsv ${genome_tsv} ..."
	}
	output {
		Map[String,String] genome = read_map(genome_tsv)
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
		if len(arr):
		    sum_ = sum(arr)
		    mean_ = sum(arr)/float(len(arr))
		    print(int(round(mean_)))
		else:
		    print(0)
		CODE
	>>>
	output {
		Int rounded_mean = read_int(stdout())
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
