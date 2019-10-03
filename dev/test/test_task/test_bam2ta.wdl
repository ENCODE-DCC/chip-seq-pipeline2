# ENCODE DCC ChIP-Seq pipeline tester for task bam2ta
# Author: Jin Lee (leepc12@gmail.com)
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bam2ta {
	Int bam2ta_subsample

	String pe_nodup_bam
	String se_nodup_bam

	String ref_pe_ta
	String ref_pe_ta_subsample
	String ref_se_ta
	String ref_se_ta_subsample
	String mito_chr_name = 'chrM'

	Int bam2ta_cpu = 1
	Int bam2ta_mem_mb = 10000
	Int bam2ta_time_hr = 6
	String bam2ta_disks = 'local-disk 100 HDD'

	call chip.bam2ta as pe_bam2ta { input :
		bam = pe_nodup_bam,
		subsample = 0,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		cpu = bam2ta_cpu,
		mem_mb = bam2ta_mem_mb,
		time_hr = bam2ta_time_hr,
		disks = bam2ta_disks,
	}
	call chip.bam2ta as pe_bam2ta_subsample { input :
		bam = pe_nodup_bam,
		subsample = bam2ta_subsample,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		cpu = bam2ta_cpu,
		mem_mb = bam2ta_mem_mb,
		time_hr = bam2ta_time_hr,
		disks = bam2ta_disks,
	}
	call chip.bam2ta as se_bam2ta { input :
		bam = se_nodup_bam,
		subsample = 0,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		cpu = bam2ta_cpu,
		mem_mb = bam2ta_mem_mb,
		time_hr = bam2ta_time_hr,
		disks = bam2ta_disks,
	}
	call chip.bam2ta as se_bam2ta_subsample { input :
		bam = se_nodup_bam,
		subsample = bam2ta_subsample,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		cpu = bam2ta_cpu,
		mem_mb = bam2ta_mem_mb,
		time_hr = bam2ta_time_hr,
		disks = bam2ta_disks,
	}

	call compare_md5sum.compare_md5sum { input :
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
