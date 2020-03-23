# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_choose_ctl {	
	String se_ta_rep1
	String se_ta_rep2
	String se_ta_pooled
	String se_ctl_ta_rep1
	String se_ctl_ta_rep2
	String se_ctl_ta_pooled

	String ref_se_choose_ctl_idx1
	String ref_se_choose_ctl_idx2
	String ref_se_choose_ctl_always_use_pooled_ctl_idx1
	String ref_se_choose_ctl_always_use_pooled_ctl_idx2
	String ref_se_choose_ctl_single_rep_idx1
	String ref_se_choose_ctl_single_ctl_idx1
	String ref_se_choose_ctl_single_ctl_idx2

	String? null

	Float ctl_depth_ratio

	Array[File] ctl_tas = [se_ctl_ta_rep1, se_ctl_ta_rep2]
	call chip.choose_ctl as se_choose_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = ctl_tas,
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	Int ctl_id_se_choose_ctl_idx1 = se_choose_ctl.chosen_ctl_ta_ids[0]
	Int ctl_id_se_choose_ctl_idx2 = se_choose_ctl.chosen_ctl_ta_ids[1]
	File ctl_se_choose_ctl_idx1 = if ctl_id_se_choose_ctl_idx1>=0 then ctl_tas[ctl_id_se_choose_ctl_idx1] else se_ctl_ta_pooled
	File ctl_se_choose_ctl_idx2 = if ctl_id_se_choose_ctl_idx2>=0 then ctl_tas[ctl_id_se_choose_ctl_idx2] else se_ctl_ta_pooled

	Array[File] ctl_tas_always_use_pooled_ctl = [se_ctl_ta_rep1, se_ctl_ta_rep2]
	call chip.choose_ctl as se_choose_ctl_always_use_pooled_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = ctl_tas_always_use_pooled_ctl,
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = true,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	Int ctl_id_se_choose_ctl_always_use_pooled_ctl_idx1 = se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_ids[0]
	Int ctl_id_se_choose_ctl_always_use_pooled_ctl_idx2 = se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_ids[1]
	File ctl_se_choose_ctl_always_use_pooled_ctl_idx1 = if ctl_id_se_choose_ctl_always_use_pooled_ctl_idx1>=0 then ctl_tas_always_use_pooled_ctl[ctl_id_se_choose_ctl_always_use_pooled_ctl_idx1] else se_ctl_ta_pooled
	File ctl_se_choose_ctl_always_use_pooled_ctl_idx2 = if ctl_id_se_choose_ctl_always_use_pooled_ctl_idx2>=0 then ctl_tas_always_use_pooled_ctl[ctl_id_se_choose_ctl_always_use_pooled_ctl_idx2] else se_ctl_ta_pooled

	Array[File] ctl_tas_single_rep = [se_ctl_ta_rep2, se_ctl_ta_rep1]
	call chip.choose_ctl as se_choose_ctl_single_rep { input :
		tas = [se_ta_rep1],
		ctl_tas = ctl_tas_single_rep,
		ta_pooled = null,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	Int ctl_id_se_choose_ctl_single_rep_idx1 = se_choose_ctl_single_rep.chosen_ctl_ta_ids[0]
	File ctl_se_choose_ctl_single_rep_idx1 = if ctl_id_se_choose_ctl_single_rep_idx1>=0 then ctl_tas_single_rep[ctl_id_se_choose_ctl_single_rep_idx1] else se_ctl_ta_pooled

	Array[File] ctl_tas_single_ctl = [se_ctl_ta_rep1]
	call chip.choose_ctl as se_choose_ctl_single_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = ctl_tas_single_ctl,
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = null,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	Int ctl_id_se_choose_ctl_single_ctl_idx1 = se_choose_ctl_single_ctl.chosen_ctl_ta_ids[0]
	Int ctl_id_se_choose_ctl_single_ctl_idx2 = se_choose_ctl_single_ctl.chosen_ctl_ta_ids[1]
	File ctl_se_choose_ctl_single_ctl_idx1 = if ctl_id_se_choose_ctl_single_ctl_idx1>=0 then ctl_tas_single_ctl[ctl_id_se_choose_ctl_single_ctl_idx1] else se_ctl_ta_pooled
	File ctl_se_choose_ctl_single_ctl_idx2 = if ctl_id_se_choose_ctl_single_ctl_idx2>=0 then ctl_tas_single_ctl[ctl_id_se_choose_ctl_single_ctl_idx2] else se_ctl_ta_pooled

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'se_choose_ctl_idx1',
			'se_choose_ctl_idx2',
			'se_choose_ctl_always_use_pooled_ctl_idx1',
			'se_choose_ctl_always_use_pooled_ctl_idx2',
			'se_choose_ctl_single_rep_idx1',
			'se_choose_ctl_single_ctl_idx1',
			'se_choose_ctl_single_ctl_idx2',
		],
		files = [
			ctl_se_choose_ctl_idx1,
			ctl_se_choose_ctl_idx2,
			ctl_se_choose_ctl_always_use_pooled_ctl_idx1,
			ctl_se_choose_ctl_always_use_pooled_ctl_idx2,
			ctl_se_choose_ctl_single_rep_idx1,
			ctl_se_choose_ctl_single_ctl_idx1,
			ctl_se_choose_ctl_single_ctl_idx2,
		],
		ref_files = [
			ref_se_choose_ctl_idx1,
			ref_se_choose_ctl_idx2,
			ref_se_choose_ctl_always_use_pooled_ctl_idx1,
			ref_se_choose_ctl_always_use_pooled_ctl_idx2,
			ref_se_choose_ctl_single_rep_idx1,
			ref_se_choose_ctl_single_ctl_idx1,
			ref_se_choose_ctl_single_ctl_idx2,
		],
	}
}

task write_to_file {
	Int i
	command {
		echo ${i} > 'tmp.txt'
	}
	output {
		File f = glob('tmp.txt')[0]
	}
}