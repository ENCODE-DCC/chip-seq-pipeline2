# ENCODE DCC ChIP-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../chipseq.wdl" as chipseq

workflow test_reproducibility {
	String se_overlap_peak_rep1_vs_rep2
	String se_overlap_peak_rep1_pr
	String se_overlap_peak_rep2_pr
	String se_overlap_peak_ppr

	String ref_se_reproducibility_qc

	call chipseq.reproducibility as se_reproducibility { input :
		prefix = 'overlap',
		peaks = [se_overlap_peak_rep1_vs_rep2],
		peaks_pr = [se_overlap_peak_rep1_pr, se_overlap_peak_rep2_pr],
		peak_ppr = se_overlap_peak_ppr,
	}

	call chipseq.compare_md5sum { input :
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
