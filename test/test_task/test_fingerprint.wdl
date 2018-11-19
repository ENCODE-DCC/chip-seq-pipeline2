# ENCODE DCC ChIP-Seq pipeline tester for task fingerprint
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_fingerprint {
	Array[File] se_nodup_bams
	File se_ctl_nodup_bam
	File se_blacklist
	Array[File] ref_se_fingerprint_logs

	Int fingerprint_cpu = 1
	Int fingerprint_mem_mb = 12000
	Int fingerprint_time_hr = 6
	String fingerprint_disks = "local-disk 100 HDD"

	call chip.fingerprint as se_fingerprint { input :
		nodup_bams = se_nodup_bams,
		ctl_bam = se_ctl_nodup_bam, # use first control only
		blacklist = se_blacklist,

		cpu = fingerprint_cpu,
		mem_mb = fingerprint_mem_mb,
		time_hr = fingerprint_time_hr,
		disks = fingerprint_disks,
	}

	# take first 8 columns (vaule in other columns are random)
	scatter(i in range(2)){
		call take_8_cols { input :
			f = se_fingerprint.jsd_qcs[i],
		}
		call take_8_cols as ref_take_8_cols { input :
			f = ref_se_fingerprint_logs[i],
		}
	}

	call chip.compare_md5sum { input :
		labels = [
			'se_fingerprint_rep1',
			'se_fingerprint_rep2',
		],
		files = [
			take_8_cols.out[0],
			take_8_cols.out[1],
		],
		ref_files = [
			ref_take_8_cols.out[0],
			ref_take_8_cols.out[1],
		],
	}
}

task take_8_cols {
	File f
	command {
		cut -f 1-8 ${f} > out.txt
	}
	output {
		File out = "out.txt"
	}
}