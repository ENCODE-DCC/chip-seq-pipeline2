# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip
import "compare_md5sum.wdl" as compare_md5sum

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

	call chip.choose_ctl as se_choose_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = [se_ctl_ta_rep1, se_ctl_ta_rep2],
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	call chip.choose_ctl as se_choose_ctl_always_use_pooled_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = [se_ctl_ta_rep1, se_ctl_ta_rep2],
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = true,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	call chip.choose_ctl as se_choose_ctl_single_rep { input :
		tas = [se_ta_rep1],
		ctl_tas = [se_ctl_ta_rep2, se_ctl_ta_rep1],
		ta_pooled = null,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
	}
	call chip.choose_ctl as se_choose_ctl_single_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = [se_ctl_ta_rep1],
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = null,
		ctl_depth_ratio = ctl_depth_ratio,
	}

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
			se_choose_ctl.chosen_ctl_tas[0],
			se_choose_ctl.chosen_ctl_tas[1],
			se_choose_ctl_always_use_pooled_ctl.chosen_ctl_tas[0],
			se_choose_ctl_always_use_pooled_ctl.chosen_ctl_tas[1],
			se_choose_ctl_single_rep.chosen_ctl_tas[0],
			se_choose_ctl_single_ctl.chosen_ctl_tas[0],
			se_choose_ctl_single_ctl.chosen_ctl_tas[1],
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