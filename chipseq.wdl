# ENCODE DCC TF/Histone ChIP-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow chipseq {
	###### pipeline inputs ######
	
	String pipeline_type  		# tf or histone
	String? peak_caller 		# force to use peak_caller if specifed. macs2 or spp

	# mandatory input files
	# experiment replicates
	Array[Array[Array[String]]]? fastqs 
								# [rep_id][merge_id][end_id] if starting from fastqs
								# 	after merging, it will reduce to 
								# 	[rep_id][end_id]
	Array[String]? bams 		# [rep_id] if starting from bams
	Array[String]? nodup_bams 	# [rep_id] if starting from filtered bams
	Array[String]? tas 			# [rep_id] if starting from tag-aligns

	Array[String]? peaks		# [rep_id] if starting from peaks
	Array[String]? peaks_pr1	# [rep_id] if starting from peaks
	Array[String]? peaks_pr2	# [rep_id] if starting from peaks
	File? peak_ppr1				# if starting from peaks
	File? peak_ppr2				# if starting from peaks
	File? peak_pooled			# if starting from peaks

	# control replicates
	Array[Array[Array[String]]]? ctl_fastqs 
									# [rep_id][merge_id][end_id]
									# 	after merging, it will reduce to 
									# 	[rep_id][end_id]
	Array[String]? ctl_bams 		# [rep_id] if starting from bams
	Array[String]? ctl_nodup_bams	# [rep_id] if starting from filtered bams
	Array[String]? ctl_tas 			# [rep_id] if starting from tag-aligns

	# mandatory genome param
	String genome_tsv 		# reference genome data TSV file including
							# all important genome specific data file paths
							# and parameters
	Boolean paired_end 		# endedness of sample

	# optional but important
	Boolean? align_only 	# disable downstream analysis (peak calling, ...)
							# after alignment
	Boolean? true_rep_only 	# disable all analyses for pseudo replicates
							# naive-overlap and IDR will also be disabled

	# task-specific variables but defined in workflow level (limit of WDL)
	# optional for MACS2
	Int? macs2_cap_num_peak	# cap number of raw peaks called from MACS2
	Float? pval_thresh 		# p.value threshold
	Int? macs2_mem_mb 		# resource (memory in MB)
	Int? macs2_time_hr		# resource (walltime in hour)
	String? macs2_disks 	# resource disks for cloud platforms

	# optional for SPP
	Int? spp_cap_num_peak	# cap number of raw peaks called from SPP
	Int? spp_cpu	 		# resource (cpu)
	Int? spp_mem_mb 		# resource (memory in MB)
	Int? spp_time_hr		# resource (walltime in hour)
	String? spp_disks 		# resource disks for cloud platforms

	# optional for IDR
	Float? idr_thresh		# IDR threshold

	# OTHER IMPORTANT mandatory/optional parameters are declared in a task level

	###### initialization for pipeline ######

	# temp null variable for optional File/String
	String? null

	# read genome data and paths
	call read_genome_tsv { input: genome_tsv = genome_tsv }	# For Google JES backend
	String bwa_idx_tar = read_genome_tsv.genome['bwa_idx_tar']
	String? blacklist = if read_genome_tsv.genome['blacklist']=='/dev/null' then null 
					else read_genome_tsv.genome['blacklist']
	String chrsz = read_genome_tsv.genome['chrsz']
	String gensz = read_genome_tsv.genome['gensz']

	# peak_caller
	String peak_caller_ = if defined(peak_caller) then select_first([peak_caller])
						else if pipeline_type=='tf' then 'spp'
						else if pipeline_type=='histone' then 'macs2'
						else 'macs2'
	# simplified variables for optional flags
	Boolean align_only_ = select_first([align_only, false])
	Boolean true_rep_only_ = select_first([true_rep_only, false])

	###### pipeline starts here ######

	Array[Array[Array[String]]] fastqs_ = select_first([fastqs, []])	
	if ( length(fastqs_)>0 ) {
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
		if ( paired_end ) {
			scatter(fastq_set in merge_fastq.merged_fastqs) {
				# for paired end dataset, map R1 only as SE for xcor analysis
				call trim_fastq { input :
					fastq = fastq_set[0],
				}
				call bwa as bwa_R1 { input :
					idx_tar = bwa_idx_tar,
					fastqs = [trim_fastq.trimmed_fastq],
					paired_end = false,
				}
				call filter as filter_R1 { input :
					bam = bwa_R1.bam,
					paired_end = false,
				}			
				call bam2ta as bam2ta_R1 { input :
					bam = filter_R1.nodup_bam,
					disable_tn5_shift = true,
					paired_end = false,
				}
			}
		}
	}
	Array[String] bams_ = select_first([bams, bwa.bam, []])
	if ( length(bams_)>0 ) {
		scatter(bam in bams_) {
			# filter/dedup bam
			call filter { input :
				bam = bam,
				paired_end = paired_end,
			}
		}
	}
	Array[String] nodup_bams_ = select_first([nodup_bams, filter.nodup_bam, []])
	if ( length(nodup_bams_)>0 ) {
		scatter(bam in nodup_bams_) {
			# convert bam to tagalign and subsample it if necessary
			call bam2ta { input :
				bam = bam,
				disable_tn5_shift = true,
				paired_end = paired_end,
			}
		}
	}
	Array[String] tas_ = select_first([tas, bam2ta.ta, []])
	Int tas_len = length(tas_)
	if ( tas_len>0 ) {
		scatter(i in range(tas_len)) {
			# for paired end dataset, if not starting from fastqs, use old method
			# (mapping with both ends for tag-aligns to be used for xcor)
			# subsample tagalign (non-mito) and cross-correlation analysis
			call xcor { input :
				ta = select_first([bam2ta_R1.ta, tas_])[i],
				paired_end = if defined(bam2ta_R1.ta) then false else paired_end,
			}
		}
		if ( !align_only_ && tas_len>1 )  {		
			# pool tagaligns from true replicates
			call pool_ta { input :
				tas = tas_,
			}
		}
	}
	# align controls
	Array[Array[Array[String]]] ctl_fastqs_ = select_first([ctl_fastqs, []])
	Int ctl_fastqs_len = length(ctl_fastqs_)
	if ( ctl_fastqs_len>0 ) {
		scatter(i in range(ctl_fastqs_len)) {
			# merge fastqs
			call merge_fastq as merge_fastq_ctl { input :
				fastqs = ctl_fastqs_[i],
				paired_end = paired_end,
			}
			# align merged fastqs with bwa
			call bwa as bwa_ctl { input :
				idx_tar = bwa_idx_tar,
				fastqs = merge_fastq_ctl.merged_fastqs, #[R1,R2]
				paired_end = paired_end,
			}
		}
	}
	Array[String] ctl_bams_ = select_first([ctl_bams, bwa_ctl.bam, []])
	if ( length(ctl_bams_)>0 ) {
		scatter(bam in ctl_bams_) {
			# filter/dedup bam
			call filter as filter_ctl { input :
				bam = bam,
				paired_end = paired_end,
			}
		}
	}
	Array[String] ctl_nodup_bams_ = select_first([ctl_nodup_bams, filter_ctl.nodup_bam, []])
	if ( length(ctl_nodup_bams_)>0 ) {
		scatter(bam in ctl_nodup_bams_) {	
			# convert bam to tagalign and subsample it if necessary
			call bam2ta as bam2ta_ctl { input :
				bam = bam,
				disable_tn5_shift = true,
				paired_end = paired_end,
			}
		}
	}
	Array[String] ctl_tas_ = select_first([ctl_tas,bam2ta_ctl.ta,[]])
	Int ctl_tas_len = length(ctl_tas_)
	if ( !align_only_ && ctl_tas_len>1 ) {	
		# pool tagaligns from true replicates
		call pool_ta as pool_ta_ctl { input :
			tas = ctl_tas_,
		}
	}

	# fingerprint and jsd plot
	if ( length(nodup_bams_)>0 && length(ctl_nodup_bams_)>0 ) {
		call fingerprint { input :
			nodup_bams = nodup_bams_,
			ctl_bam = ctl_nodup_bams_[0], # use first control only
			blacklist = blacklist,
		}
	}

	# call peaks because we have all tas and ctl_tas (optional for histone chipseq) ready
	if ( !align_only_ && tas_len>0 ) {
		# choose appropriate control for each exp IP replicate
		# outputs:
		# 	choose_ctl.idx : control replicate index for each exp replicate 
		#					-1 means pooled ctl replicate
		call choose_ctl { input:
			tas = tas_,
			ctl_tas = ctl_tas_,
			ta_pooled = if tas_len>1 then pool_ta.ta_pooled else null,
			ctl_ta_pooled = if ctl_tas_len>1 then pool_ta_ctl.ta_pooled else null,
		}
		scatter(i in range(tas_len)) {
			# choose appropriate control based on index (choose_ctl.idx)
			String? ctl_ta = if ctl_tas_len==0 then null
					else if choose_ctl.idx[i]<0 then pool_ta_ctl.ta_pooled
					else ctl_tas_[(choose_ctl.idx[i])]
		}
		scatter(i in range(tas_len)) {
			# always call MACS2 peaks for true replicates to get signal tracks
			# call peaks on tagalign
			call macs2 { input :
				ta = tas_[i],
				ctl_ta = ctl_ta[i],
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				make_signal = true,
				fraglen = select_first([xcor.fraglen])[i],
				blacklist = blacklist,
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
		if ( peak_caller_=='spp' ) {
			scatter(i in range(tas_len)) {
				# call peaks on tagalign
				call spp { input :
					ta = tas_[i],
					ctl_ta = ctl_ta[i],
					chrsz = chrsz,
					cap_num_peak = spp_cap_num_peak,
					fraglen = select_first([xcor.fraglen])[i],
					blacklist = blacklist,
					cpu = spp_cpu,
					mem_mb = spp_mem_mb,
					disks = spp_disks,
					time_hr = spp_time_hr,
				}
			}
		}
		if ( !true_rep_only_ ) {
			scatter(i in range(tas_len)) {
				# make two self pseudo replicates per true replicate
				call spr { input :
					ta = tas_[i],
					paired_end = paired_end,
				}				
			}			
			if ( peak_caller_=='macs2' ) {
				scatter(i in range(tas_len)) {				
					# call peaks on 1st pseudo replicated tagalign 
					call macs2 as macs2_pr1 { input :
						ta = spr.ta_pr1[i],
						ctl_ta = ctl_ta[i],
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = macs2_cap_num_peak,
						pval_thresh = pval_thresh,
						fraglen = select_first([xcor.fraglen])[i],
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}
					# call peaks on 2nd pseudo replicated tagalign 
					call macs2 as macs2_pr2 { input :
						ta = spr.ta_pr2[i],
						ctl_ta = ctl_ta[i],
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = macs2_cap_num_peak,
						pval_thresh = pval_thresh,
						fraglen = select_first([xcor.fraglen])[i],
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}					
				}
			}
			if ( peak_caller_=='spp' ) {
				scatter(i in range(tas_len)) {				
					# call peaks on 1st pseudo replicated tagalign 
					call spp as spp_pr1 { input :
						ta = spr.ta_pr1[i],
						ctl_ta = ctl_ta[i],
						chrsz = chrsz,
						cap_num_peak = spp_cap_num_peak,
						fraglen = select_first([xcor.fraglen])[i],
						blacklist = blacklist,
						cpu = spp_cpu,
						mem_mb = spp_mem_mb,
						disks = spp_disks,
						time_hr = spp_time_hr,
					}
					# call peaks on 2nd pseudo replicated tagalign 
					call spp as spp_pr2 { input :
						ta = spr.ta_pr2[i],
						ctl_ta = ctl_ta[i],
						chrsz = chrsz,
						cap_num_peak = spp_cap_num_peak,
						fraglen = select_first([xcor.fraglen])[i],
						blacklist = blacklist,
						cpu = spp_cpu,
						mem_mb = spp_mem_mb,
						disks = spp_disks,
						time_hr = spp_time_hr,
					}
				}
			}
		}
		if ( tas_len>1 ) {
			String? ctl_ta_pooled = if ctl_tas_len==0 then null
						else if ctl_tas_len>1 then pool_ta_ctl.ta_pooled
						else ctl_tas_[0]
			# get sum, rounded mean of fragment length
			# these will be used for 
			# 1) calling peaks for pooled true/pseudo replicates
			# 2) calculating FRiP
			call rounded_mean { input :
				ints = select_first([xcor.fraglen]),
			}
			# call peaks on pooled replicate
			# always call MACS2 peaks for pooled replicate to get signal tracks
			call macs2 as macs2_pooled { input :
				ta = pool_ta.ta_pooled,
				ctl_ta = ctl_ta_pooled,
				gensz = gensz,
				chrsz = chrsz,
				cap_num_peak = macs2_cap_num_peak,
				pval_thresh = pval_thresh,
				make_signal = true,
				fraglen = rounded_mean.rounded_mean,
				blacklist = blacklist,
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
			if ( peak_caller_=='spp' ) {
				# call peaks on pooled replicate
				call spp as spp_pooled { input :
					ta = pool_ta.ta_pooled,
					ctl_ta = ctl_ta_pooled,
					chrsz = chrsz,
					cap_num_peak = spp_cap_num_peak,
					fraglen = rounded_mean.rounded_mean,
					blacklist = blacklist,
					cpu = spp_cpu,
					mem_mb = spp_mem_mb,
					disks = spp_disks,
					time_hr = spp_time_hr,
				}
			}
			if ( !true_rep_only_ ) {
				# pool tagaligns from pseudo replicates
				call pool_ta as pool_ta_pr1 { input :
					tas = spr.ta_pr1,
				}
				call pool_ta as pool_ta_pr2 { input :
					tas = spr.ta_pr2,
				}				
				if ( peak_caller_=='macs2' ) {
					# call peaks on 1st pooled pseudo replicates
					call macs2 as macs2_ppr1 { input :
						ta = pool_ta_pr1.ta_pooled,
						ctl_ta = ctl_ta_pooled,
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = macs2_cap_num_peak,
						pval_thresh = pval_thresh,
						fraglen = rounded_mean.rounded_mean,
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}
					# call peaks on 2nd pooled pseudo replicates
					call macs2 as macs2_ppr2 { input :
						ta = pool_ta_pr2.ta_pooled,
						ctl_ta = ctl_ta_pooled,
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = macs2_cap_num_peak,
						pval_thresh = pval_thresh,
						fraglen = rounded_mean.rounded_mean,
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}
				}
				if ( peak_caller_=='spp' ) {
					# call peaks on 1st pooled pseudo replicates
					call spp as spp_ppr1 { input :
						ta = pool_ta_pr1.ta_pooled,
						ctl_ta = ctl_ta_pooled,
						chrsz = chrsz,
						cap_num_peak = spp_cap_num_peak,
						fraglen = rounded_mean.rounded_mean,
						blacklist = blacklist,
						cpu = spp_cpu,
						mem_mb = spp_mem_mb,
						disks = spp_disks,
						time_hr = spp_time_hr,
					}
					# call peaks on 2nd pooled pseudo replicates
					call spp as spp_ppr2 { input :
						ta = pool_ta_pr2.ta_pooled,
						ctl_ta = ctl_ta_pooled,
						chrsz = chrsz,
						cap_num_peak = spp_cap_num_peak,
						fraglen = rounded_mean.rounded_mean,
						blacklist = blacklist,
						cpu = spp_cpu,
						mem_mb = spp_mem_mb,
						disks = spp_disks,
						time_hr = spp_time_hr,
					}
				}
			}
		}
	}

	Array[String?] peaks_ = select_first([peaks, spp.rpeak, macs2.npeak, []])
	Array[String?] peaks_pr1_ = select_first([peaks_pr1, spp_pr1.rpeak, macs2_pr1.npeak, []])
	Array[String?] peaks_pr2_ = select_first([peaks_pr2, spp_pr2.rpeak, macs2_pr2.npeak, []])
	Int num_rep = length(peaks_)
	String? peak_pooled_ = select_first([peak_pooled, spp_pooled.rpeak, macs2_pooled.npeak, '/dev/null'])
	String? peak_ppr1_ = select_first([peak_ppr1, spp_ppr1.rpeak, macs2_ppr1.npeak, '/dev/null'])
	String? peak_ppr2_ = select_first([peak_ppr2, spp_ppr2.rpeak, macs2_ppr2.npeak, '/dev/null'])
	# determine peak_type
	String peak_type = if peak_caller_=='macs2' then 'narrowPeak'
					else if peak_caller_=='spp' then 'regionPeak'
					else 'narrowPeak'
	# determine idr ranking method
	String idr_rank = if peak_caller_=='macs2' then 'p.value'
					else if peak_caller_=='spp' then 'signal.value'
					else 'p.value'
	# enable_idr for TF chipseq only
	Boolean enable_idr = pipeline_type=='tf'
	
	# generate all possible pairs of true replicates
	call pair_gen { input: num_rep = num_rep }

	if ( !align_only_ ) {
		if ( num_rep>1 ) {
			# Naive overlap on every pair of true replicates
			scatter( pair in pair_gen.pairs ) {
				call overlap { input :
					prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
					peak1 = peaks_[(pair[0])],
					peak2 = peaks_[(pair[1])],
					peak_pooled = peak_pooled_,
					peak_type = peak_type,
					blacklist = blacklist,
					ta = if length(tas_)>0 then pool_ta.ta_pooled else null,
				}
			}
			if ( enable_idr ) {
				# IDR on every pair of true replicates
				scatter( pair in pair_gen.pairs ) {
					call idr { input : 
						prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
						peak1 = peaks_[(pair[0])],
						peak2 = peaks_[(pair[1])],
						peak_pooled = peak_pooled_,
						idr_thresh = select_first([idr_thresh,0.1]),
						peak_type = peak_type,
						rank = idr_rank,
						blacklist = blacklist,
						ta = if length(tas_)>0 then pool_ta.ta_pooled else null,
					}
				}
			}
		}
		if ( !true_rep_only_ ) {
			# Naive overlap on pseduo replicates
			scatter( i in range(num_rep) ) {
				call overlap as overlap_pr { input : 
					prefix = "rep"+(i+1)+"-pr",
					peak1 = peaks_pr1_[i],
					peak2 = peaks_pr2_[i],
					peak_pooled = peaks_[i],
					peak_type = peak_type,
					blacklist = blacklist,
					ta = if length(tas_)>0 then tas_[i] else null,
				}
			}
			if ( enable_idr ) {
				# IDR on pseduo replicates
				scatter( i in range(num_rep) ) {
					call idr as idr_pr { input : 
						prefix = "rep"+(i+1)+"-pr",
						peak1 = peaks_pr1_[i],
						peak2 = peaks_pr2_[i],
						peak_pooled = peaks_[i],
						idr_thresh = select_first([idr_thresh,0.1]),
						peak_type = peak_type,
						rank = idr_rank,
						blacklist = blacklist,
						ta = if length(tas_)>0 then tas_[i] else null,
					}
				}
			}
			if ( num_rep>1 ) {
				# Naive overlap on pooled pseudo replicates
				call overlap as overlap_ppr { input : 
					prefix = "ppr",
					peak1 = peak_ppr1_,
					peak2 = peak_ppr2_,
					peak_pooled = peak_pooled_,
					peak_type = peak_type,
					blacklist = blacklist,
					ta = if length(tas_)>0 then pool_ta.ta_pooled else null,
				}
				if ( enable_idr ) {
					# IDR on pooled pseduo replicates
					call idr as idr_ppr { input : 
						prefix = "ppr",
						peak1 = peak_ppr1_,
						peak2 = peak_ppr2_,
						peak_pooled = peak_pooled_,
						idr_thresh = select_first([idr_thresh,0.1]),
						peak_type = peak_type,
						rank = idr_rank,
						blacklist = blacklist,
						ta = if length(tas_)>0 then pool_ta.ta_pooled else null,
					}
				}
			}
			# reproducibility QC for overlapping peaks
			call reproducibility as reproducibility_overlap { input :
				prefix = 'overlap',
				peaks = select_first([overlap.bfilt_overlap_peak, []]),
				peaks_pr = overlap_pr.bfilt_overlap_peak,
				peak_ppr = overlap_ppr.bfilt_overlap_peak,
			}
			if ( enable_idr ) {
				# reproducibility QC for IDR peaks
				call reproducibility as reproducibility_idr { input :
					prefix = 'idr',
					peaks = select_first([idr.bfilt_idr_peak, []]),
					peaks_pr = idr_pr.bfilt_idr_peak,
					peak_ppr = idr_ppr.bfilt_idr_peak,
				}
			}
		}
	}

	call qc_report { input :
		paired_end = paired_end,
		pipeline_type = pipeline_type,
		peak_caller = peak_caller_,
		idr_thresh = select_first([idr_thresh,0.1]),

		flagstat_qcs = select_first([bwa.flagstat_qc, []]),
		nodup_flagstat_qcs = select_first([filter.flagstat_qc, []]),
		dup_qcs = select_first([filter.dup_qc, []]),
		pbc_qcs = select_first([filter.pbc_qc, []]),
		
		ctl_flagstat_qcs = select_first([bwa_ctl.flagstat_qc, []]),
		ctl_nodup_flagstat_qcs = select_first([filter_ctl.flagstat_qc, []]),
		ctl_dup_qcs = select_first([filter_ctl.dup_qc, []]),
		ctl_pbc_qcs = select_first([filter_ctl.pbc_qc, []]),

		xcor_plots = select_first([xcor.plot_png, []]),
		xcor_scores = select_first([xcor.score, []]),

		jsd_plot = if defined(fingerprint.plot) then fingerprint.plot else null,
		jsd_qcs = select_first([fingerprint.jsd_qcs, []]),

		frip_qcs = select_first([spp.frip_qc, macs2.frip_qc, []]),
		frip_qcs_pr1 = select_first([spp_pr1.frip_qc, macs2_pr1.frip_qc, []]),
		frip_qcs_pr2 = select_first([spp_pr2.frip_qc, macs2_pr2.frip_qc, []]),
		frip_qc_pooled = if defined(spp_pooled.frip_qc) then spp_pooled.frip_qc
						else if defined(macs2_pooled.frip_qc) then macs2_pooled.frip_qc
						else null,
		frip_qc_ppr1 = if defined(spp_ppr1.frip_qc) then spp_ppr1.frip_qc
						else if defined(macs2_ppr1.frip_qc) then macs2_ppr1.frip_qc
						else null,
		frip_qc_ppr2 = if defined(spp_ppr2.frip_qc) then spp_ppr2.frip_qc
						else if defined(macs2_ppr2.frip_qc) then macs2_ppr2.frip_qc
						else null,

		idr_plots = select_first([idr.idr_plot, []]),
		idr_plots_pr = select_first([idr_pr.idr_plot, []]),
		idr_plot_ppr = if defined(idr_ppr.idr_plot) then idr_ppr.idr_plot else null,
		frip_idr_qcs = select_first([idr.frip_qc, []]),
		frip_idr_qcs_pr = select_first([idr_pr.frip_qc, []]),
		frip_idr_qc_ppr = if defined(idr_ppr.frip_qc) then idr_ppr.frip_qc else null,
		frip_overlap_qcs = select_first([overlap.frip_qc, []]),
		frip_overlap_qcs_pr = select_first([overlap_pr.frip_qc, []]),
		frip_overlap_qc_ppr = if defined(overlap_ppr.frip_qc) then overlap_ppr.frip_qc else null,
		idr_reproducibility_qc = if defined(reproducibility_idr.reproducibility_qc) 
								then reproducibility_idr.reproducibility_qc else null,
		overlap_reproducibility_qc = if defined(reproducibility_overlap.reproducibility_qc) 
								then reproducibility_overlap.reproducibility_qc else null,
	}
}

### genomic tasks

task merge_fastq { # trim adapters and merge trimmed fastqs
	# parameters from workflow
	Array[Array[File]] fastqs 		# [merge_id][end_id]
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
}

task bwa {
	# parameters from workflow
	File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [end_id]
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

task fingerprint {
	# parameters from workflow
	Array[File?] nodup_bams
	File? ctl_bam	 		# one control bam is required
	File? blacklist

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

task xcor {
	# parameters from workflow
	File ta
	Boolean paired_end
	# optional
	Int? subsample 		# number of reads to subsample TAGALIGN
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
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task choose_ctl {
	# parameters from workflow
	Array[File?] tas
	Array[File?] ctl_tas
	File? ta_pooled
	File? ctl_ta_pooled
	# optional
	Boolean? always_use_pooled_ctl # always use pooled control for all exp rep.
	Float? ctl_depth_ratio 	# if ratio between controls is higher than this
							# then always use pooled control for all exp rep.
	command {
		python $(which encode_choose_ctl.py) \
			--tas ${sep=' ' tas} \
			--ctl-tas ${sep=' ' ctl_tas} \
			${"--ta-pooled " + ta_pooled} \
			${"--ctl-ta-pooled " + ctl_ta_pooled} \
			${if select_first([always_use_pooled_ctl,false]) then 
				"--always-use-pooled-ctl" else ""} \
			${"--ctl-depth-ratio " + select_first([ctl_depth_ratio,"1.2"])}				
	}
	output {
		Array[Int] idx = read_lines("idx.txt")
	}
}

task macs2 {
	# parameters from workflow
	File? ta
	File? ctl_ta 		# optional. macs2 can work without control
	Int? fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	Int? cap_num_peak	# cap number of raw peaks called from MACS2
	Float? pval_thresh	# p.value threshold
	Boolean? make_signal
	File? blacklist 	# blacklist BED to filter raw peaks
	# fixed var
	String peak_type = "narrowPeak"
	# resource
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_macs2_chipseq.py) \
			${ta} \
			${"--ctl-ta " + ctl_ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--cap-num-peak " + select_first([cap_num_peak,500000])} \
			${"--pval-thresh "+ pval_thresh} \
			${if select_first([make_signal,false]) then "--make-signal" else ""} \
			${"--blacklist "+ blacklist}

		${if select_first([make_signal,false]) then "" 
			else "touch null.pval.signal.bigwig null.fc.signal.bigwig"}
		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}			
		touch null # ugly part to deal with optional outputs
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_npeak = if defined(blacklist) then glob("*.bfilt."+peak_type+".gz")[0] else npeak
		File sig_pval = if select_first([make_signal,false]) then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if select_first([make_signal,false]) then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task spp {
	# parameters from workflow
	File? ta
	File? ctl_ta
	Int? fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	Int? cap_num_peak	# cap number of raw peaks called from MACS2
	File? blacklist 	# blacklist BED to filter raw peaks
	# fixed var
	String peak_type = "regionPeak"
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_spp.py) \
			${ta} \
			${"--ctl-ta " + ctl_ta} \
			${"--chrsz " + chrsz} \
			${"--fraglen " + fraglen} \
			${"--cap-num-peak " + select_first([cap_num_peak,300000])} \
			${"--nth " + select_first([cpu,2])} \
			${"--blacklist "+ blacklist}

		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}
	}
	output {
		File rpeak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_rpeak = if defined(blacklist) then glob("*.bfilt."+peak_type+".gz")[0] else rpeak
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

task filter {
	# parameters from workflow
	File bam
	Boolean paired_end
	# optional
	String? dup_marker 			# picard.jar MarkDuplicates (picard) or 
								# sambamba markdup (sambamba)
	Int? mapq_thresh			# threshold for low MAPQ reads removal
	Boolean? no_dup_removal 	# no dupe reads removal when filtering BAM
								# dup.qc and pbc.qc will be emptry files
								# and nodup_bam in the output is 
	# resource					# filtered bam with dupes	
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if select_first([no_dup_removal,false]) then "--no-dup-removal" else ""} \
			${"--nth " + cpu}
		# ugly part to deal with optional outputs with Google JES backend
		${if select_first([no_dup_removal,false]) then 
			"touch null.dup.qc null.pbc.qc; " else ""}
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
	String? regex_grep_v_ta 	# Perl-style regular expression pattern 
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
			${"--regex-grep-v-ta " +"'"+regex_grep_v_ta+"'"} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}
	runtime {
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
		memory : "${select_first([mem_mb,'12000'])} MB"
	}
}

task pool_ta {
	# parameters from workflow
	Array[File]? tas

	command {
		python $(which encode_pool_ta.py) \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
	}
}

task idr {
	# parameters from workflow
	String? prefix 		# prefix for IDR output file
	File? peak1 			
	File? peak2
	File? peak_pooled
	Float? idr_thresh
	File? blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int? fraglen 		# fragment length from xcor
	File? chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}
		${if defined(ta) then "" 
			else "touch null.frip.qc"}			
		touch null 
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = if defined(blacklist) then 
							glob("*.bfilt."+peak_type+".gz")[0] else idr_peak
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
}

task overlap {
	# parameters from workflow
	String prefix 		# prefix for IDR output file
	File? peak1
	File? peak2
	File? peak_pooled
	File? blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int? fraglen 		# fragment length from xcor
	File? chrsz			# 2-col chromosome sizes file
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
		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}
		${if defined(ta) then "" 
			else "touch null.frip.qc"}			
		touch null 
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = if defined(blacklist) then 
							glob("*.bfilt."+peak_type+".gz")[0] else overlap_peak
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
}

task reproducibility {
	# parameters from workflow
	String prefix
	Array[File] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File]? peaks_pr	# peak files from pseudo replicates
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
	Array[File?] idr_plots
	Array[File?] idr_plots_pr
	File? idr_plot_ppr
	Array[File?] frip_qcs
	Array[File?] frip_qcs_pr1
	Array[File?] frip_qcs_pr2
	File? frip_qc_pooled
	File? frip_qc_ppr1
	File? frip_qc_ppr2
	Array[File?] frip_idr_qcs
	Array[File?] frip_idr_qcs_pr
	File? frip_idr_qc_ppr
	Array[File?] frip_overlap_qcs
	Array[File?] frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

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
			--frip-qcs ${sep=' ' frip_qcs} \
			--frip-qcs-pr1 ${sep=' ' frip_qcs_pr1} \
			--frip-qcs-pr2 ${sep=' ' frip_qcs_pr2} \
			${"--frip-qc-pooled " + frip_qc_pooled} \
			${"--frip-qc-ppr1 " + frip_qc_ppr1} \
			${"--frip-qc-ppr2 " + frip_qc_ppr2} \
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
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		String qc_json_str = read_string(qc_json)
		#File encode_accession_json= glob('*encode_accession.json')[0]
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
}

task pair_gen {
	Int num_rep
	command <<<
		python <<CODE
		for i in range(${num_rep}):
		    for j in range(i+1,${num_rep}):
		        print('{}\t{}'.format(i,j))
		CODE
	>>>
	output {
		Array[Array[Int]] pairs = if num_rep>1 then read_tsv(stdout()) else [[]]
	}
}

task rounded_mean {
	Array[Int?] ints
	command <<<
		python <<CODE
		arr = [${sep=',' ints}]
		sum_ = sum(arr)
		mean_ = sum(arr)/float(len(arr))
		print(int(round(mean_)))
		CODE
	>>>
	output {
		Int rounded_mean = read_int(stdout())
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
}
