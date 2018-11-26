# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_reproducibility {
	String se_overlap_peak_rep1_vs_rep2
	String se_overlap_peak_rep1_pr
	String se_overlap_peak_rep2_pr
	String se_overlap_peak_ppr
	String se_chrsz

	String ref_se_reproducibility_qc

	call chip.reproducibility as se_reproducibility { input :
		prefix = 'overlap',
		peaks = [se_overlap_peak_rep1_vs_rep2],
		peaks_pr = [se_overlap_peak_rep1_pr, se_overlap_peak_rep2_pr],
		peak_ppr = se_overlap_peak_ppr,
		peak_type = 'regionPeak',		
		chrsz = se_chrsz,		
		keep_irregular_chr_in_bfilt_peak = false,
	}

	call chip.compare_md5sum { input :
		labels = [
			'se_reproducibility',
		],
		files = [
			se_reproducibility.reproducibility_qc,
		],
		ref_files = [
			ref_se_reproducibility_qc,
		],
	}
}
