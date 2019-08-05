# ENCODE DCC ChIP-Seq pipeline tester for task macs2
# Author: Jin Lee (leepc12@gmail.com)
import "../../../chip.wdl" as chip
import "compare_md5sum.wdl" as compare_md5sum

workflow test_macs2_signal_track {
	Float pval_thresh

	Int fraglen
	# test macs2 for SE set only
	String se_ta
	String se_ctl_ta

	String ref_se_macs2_pval_bw # p-val signal
	String se_chrsz
	String se_gensz

	Int macs2_mem_mb = 16000
	Int macs2_time_hr = 24
	String macs2_disks = "local-disk 100 HDD"	

	call chip.macs2_signal_track as se_macs2_signal_track { input :
		tas = [se_ta, se_ctl_ta],
		gensz = se_gensz,
		chrsz = se_chrsz,
		fraglen = fraglen,
		pval_thresh = pval_thresh,

		mem_mb = macs2_mem_mb,
		time_hr = macs2_time_hr,
		disks = macs2_disks,		
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'se_macs2_pval_bw',
		],
		files = [
			se_macs2_signal_track.pval_bw,
		],
		ref_files = [
			ref_se_macs2_pval_bw,
		],
	}
}
