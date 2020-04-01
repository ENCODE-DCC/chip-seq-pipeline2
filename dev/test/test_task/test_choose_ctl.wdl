# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_choose_ctl {	
	# See ./test_choose_ctl.aux.xlsx for details
	String se_ta_rep1
	String se_ta_rep2
	String se_ta_pooled
	String se_ctl_ta_rep1
	String se_ctl_ta_rep2
	String se_ctl_ta_pooled

	Int ref_se_choose_ctl_idx1
	Int ref_se_choose_ctl_idx2
	Int ref_se_choose_ctl_sub1
	Int ref_se_choose_ctl_sub2
	Int ref_se_choose_ctl_sub_pooled

	Int ref_se_choose_ctl_always_use_pooled_ctl_idx1
	Int ref_se_choose_ctl_always_use_pooled_ctl_idx2
	Int ref_se_choose_ctl_always_use_pooled_ctl_sub1
	Int ref_se_choose_ctl_always_use_pooled_ctl_sub2
	Int ref_se_choose_ctl_always_use_pooled_ctl_sub_pooled

	Int ref_se_choose_ctl_single_rep_idx1
	Int ref_se_choose_ctl_single_rep_sub1
	Int ref_se_choose_ctl_single_rep_sub_pooled

	Int ref_se_choose_ctl_single_ctl_idx1
	Int ref_se_choose_ctl_single_ctl_idx2
	Int ref_se_choose_ctl_single_ctl_sub1
	Int ref_se_choose_ctl_single_ctl_sub2
	Int ref_se_choose_ctl_single_ctl_sub_pooled

	String? null

	Float ctl_depth_ratio
	Int ctl_depth_limit
	Float exp_ctl_depth_ratio_limit

	call chip.choose_ctl as se_choose_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = [se_ctl_ta_rep1, se_ctl_ta_rep2],
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
		ctl_depth_limit = ctl_depth_limit,
		exp_ctl_depth_ratio_limit = exp_ctl_depth_ratio_limit,
	}

	call chip.choose_ctl as se_choose_ctl_always_use_pooled_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = [se_ctl_ta_rep1, se_ctl_ta_rep2],
		ta_pooled = se_ta_pooled,
		always_use_pooled_ctl = true,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
		ctl_depth_limit = ctl_depth_limit,
		exp_ctl_depth_ratio_limit = exp_ctl_depth_ratio_limit,
	}

	# INTENTIONALLY swapped ctl_tas 
	call chip.choose_ctl as se_choose_ctl_single_rep { input :
		tas = [se_ta_rep1],
		ctl_tas = [se_ctl_ta_rep2, se_ctl_ta_rep1],
		ta_pooled = null,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = se_ctl_ta_pooled,
		ctl_depth_ratio = ctl_depth_ratio,
		ctl_depth_limit = ctl_depth_limit,
		exp_ctl_depth_ratio_limit = exp_ctl_depth_ratio_limit,
	}

	call chip.choose_ctl as se_choose_ctl_single_ctl { input :
		tas = [se_ta_rep1, se_ta_rep2],
		ctl_tas = [se_ctl_ta_rep1],
		ta_pooled = se_ctl_ta_rep1,
		always_use_pooled_ctl = false,
		ctl_ta_pooled = null,
		ctl_depth_ratio = ctl_depth_ratio,
		ctl_depth_limit = ctl_depth_limit,
		exp_ctl_depth_ratio_limit = exp_ctl_depth_ratio_limit,
	}

	Map[String, Pair[Int, Int]] tests = {
		'se_choose_ctl_idx1': (se_choose_ctl.chosen_ctl_ta_ids[0], ref_se_choose_ctl_idx1),
		'se_choose_ctl_idx2': (se_choose_ctl.chosen_ctl_ta_ids[1], ref_se_choose_ctl_idx2),
		'se_choose_ctl_sub1': (se_choose_ctl.chosen_ctl_ta_subsample[0], ref_se_choose_ctl_sub1),
		'se_choose_ctl_sub2': (se_choose_ctl.chosen_ctl_ta_subsample[1], ref_se_choose_ctl_sub2),
		'se_choose_ctl_sub_pooled': (se_choose_ctl.chosen_ctl_ta_subsample_pooled, ref_se_choose_ctl_sub_pooled),

		'se_choose_ctl_always_use_pooled_ctl_idx1': (se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_ids[0], ref_se_choose_ctl_always_use_pooled_ctl_idx1),
		'se_choose_ctl_always_use_pooled_ctl_idx2': (se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_ids[1], ref_se_choose_ctl_always_use_pooled_ctl_idx2),
		'se_choose_ctl_always_use_pooled_ctl_sub1': (se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_subsample[0], ref_se_choose_ctl_always_use_pooled_ctl_sub1),
		'se_choose_ctl_always_use_pooled_ctl_sub2': (se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_subsample[1], ref_se_choose_ctl_always_use_pooled_ctl_sub2),
		'se_choose_ctl_always_use_pooled_ctl_sub_pooled': (se_choose_ctl_always_use_pooled_ctl.chosen_ctl_ta_subsample_pooled, ref_se_choose_ctl_always_use_pooled_ctl_sub_pooled),

		'se_choose_ctl_single_rep_idx1': (se_choose_ctl_single_rep.chosen_ctl_ta_ids[0], ref_se_choose_ctl_single_rep_idx1),
		'se_choose_ctl_single_rep_sub1': (se_choose_ctl_single_rep.chosen_ctl_ta_subsample[0], ref_se_choose_ctl_single_rep_sub1),
		'se_choose_ctl_single_rep_sub_pooled': (se_choose_ctl_single_rep.chosen_ctl_ta_subsample_pooled, ref_se_choose_ctl_single_rep_sub_pooled),

		'se_choose_ctl_single_ctl_idx1': (se_choose_ctl_single_ctl.chosen_ctl_ta_ids[0], ref_se_choose_ctl_single_ctl_idx1),
		'se_choose_ctl_single_ctl_idx2': (se_choose_ctl_single_ctl.chosen_ctl_ta_ids[1], ref_se_choose_ctl_single_ctl_idx2),
		'se_choose_ctl_single_ctl_sub1': (se_choose_ctl_single_ctl.chosen_ctl_ta_subsample[0], ref_se_choose_ctl_single_ctl_sub1),
		'se_choose_ctl_single_ctl_sub2': (se_choose_ctl_single_ctl.chosen_ctl_ta_subsample[1], ref_se_choose_ctl_single_ctl_sub2),
		'se_choose_ctl_single_ctl_sub_pooled': (se_choose_ctl_single_ctl.chosen_ctl_ta_subsample_pooled, ref_se_choose_ctl_single_ctl_sub_pooled),
	}
	scatter( test in tests ) {
		String k = test.left
		Pair[Int, Int] v = test.right
		if ( v.left != v.right ) {
			call chip.raise_exception { input:
				msg = k,
				vals = [v.left, v.right]
			}
		}
	}
}
