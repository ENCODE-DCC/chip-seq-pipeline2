# ENCODE DCC ChIP-Seq pipeline tester for task filter
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_filter {
	String pe_bam
	String se_bam

	String ref_pe_nodup_bam
	String ref_pe_filt_bam
	String ref_se_nodup_bam
	String ref_se_filt_bam

	call chip.filter as pe_filter { input :
		bam = pe_bam,
		paired_end = true,
		cpu = 1,
	}
	call chip.filter as pe_filter_no_dup_removal { input :
		bam = pe_bam,
		no_dup_removal = true,
		paired_end = true,
		cpu = 1,
	}
	call chip.filter as se_filter { input :
		bam = se_bam,
		paired_end = false,
		cpu = 1,
	}
	call chip.filter as se_filter_no_dup_removal { input :
		bam = se_bam,
		no_dup_removal = true,
		paired_end = false,
		cpu = 1,
	}

	call chip.compare_md5sum { input :
		labels = [
			'pe_filter',
			'pe_filter_no_dup_removal',
			'se_filter',
			'se_filter_no_dup_removal',
		],
		files = [
			pe_filter.nodup_bam,
			pe_filter_no_dup_removal.nodup_bam,
			se_filter.nodup_bam,
			se_filter_no_dup_removal.nodup_bam,
		],
		ref_files = [
			ref_pe_nodup_bam,
			ref_pe_filt_bam,
			ref_se_nodup_bam,
			ref_se_filt_bam,
		],
	}
}
