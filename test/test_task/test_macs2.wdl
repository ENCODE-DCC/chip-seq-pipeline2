# ENCODE DCC ChIP-Seq pipeline tester for task macs2
# Author: Jin Lee (leepc12@gmail.com)
import "../../chipseq.wdl" as chipseq

workflow test_macs2 {
	Int cap_num_peak
	Float pval_thresh

	Int fraglen
	# test macs2 for SE set only
	String se_ta
	String se_ctl_ta

	String ref_se_macs2_npeak # raw narrow-peak
	String ref_se_macs2_bfilt_npeak # blacklist filtered narrow-peak
	String ref_se_macs2_frip_qc 
	String ref_se_macs2_sig_pval # p-val signal

	String se_blacklist
	String se_chrsz
	String se_gensz

	call chipseq.macs2 as se_macs2 { input :
		ta = se_ta,
		ctl_ta = se_ctl_ta,
		gensz = se_gensz,
		chrsz = se_chrsz,
		fraglen = fraglen,
		cap_num_peak = cap_num_peak,
		pval_thresh = pval_thresh,
		make_signal = true,
		blacklist = se_blacklist,
	}

	call chipseq.compare_md5sum { input :
		labels = [
			'se_macs2_npeak',
			'se_macs2_bfilt_npeak',
			'se_macs2_frip_qc',
			'se_macs2_sig_pval',
		],
		files = [
			se_macs2.npeak,
			se_macs2.bfilt_npeak,
			se_macs2.frip_qc,
			se_macs2.sig_pval,
		],
		ref_files = [
			ref_se_macs2_npeak,
			ref_se_macs2_bfilt_npeak,
			ref_se_macs2_frip_qc,
			ref_se_macs2_sig_pval,
		],
	}
}
