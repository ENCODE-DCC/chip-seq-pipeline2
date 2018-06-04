# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_pool_ta {
	String se_ta_rep1
	String se_ta_rep2

	String ref_se_pooled_ta

	call chip.pool_ta as se_pool_ta { input :
		tas = [se_ta_rep1, se_ta_rep2],
	}

	call chip.compare_md5sum { input :
		labels = [
			'se_pool_ta',
		],
		files = [
			se_pool_ta.ta_pooled,
		],
		ref_files = [
			ref_se_pooled_ta,
		],
	}
}
