# ENCODE DCC ChIP-Seq pipeline tester for task bam2ta
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_bam2ta {
	Int bam2ta_subsample

	String pe_nodup_bam
	String se_nodup_bam

	String ref_pe_ta
	String ref_pe_ta_subsample
	String ref_se_ta
	String ref_se_ta_subsample

	call chip.bam2ta as pe_bam2ta { input :
		bam = pe_nodup_bam,
		disable_tn5_shift = true,
		paired_end = true,
	}
	call chip.bam2ta as pe_bam2ta_subsample { input :
		bam = pe_nodup_bam,
		disable_tn5_shift = true,
		subsample = bam2ta_subsample,
		paired_end = true,
	}
	call chip.bam2ta as se_bam2ta { input :
		bam = se_nodup_bam,
		disable_tn5_shift = true,
		paired_end = false,
	}
	call chip.bam2ta as se_bam2ta_subsample { input :
		bam = se_nodup_bam,
		disable_tn5_shift = true,
		subsample = bam2ta_subsample,
		paired_end = false,
	}

	call chip.compare_md5sum { input :
		labels = [
			'pe_bam2ta',
			'pe_bam2ta_subsample',
			'se_bam2ta',
			'se_bam2ta_subsample',
		],
		files = [
			pe_bam2ta.ta,
			pe_bam2ta_subsample.ta,
			se_bam2ta.ta,
			se_bam2ta_subsample.ta,
		],
		ref_files = [
			ref_pe_ta,
			ref_pe_ta_subsample,
			ref_se_ta,
			ref_se_ta_subsample,
		],
	}
}
