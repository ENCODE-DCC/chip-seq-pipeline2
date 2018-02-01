# ENCODE DCC ChIP-Seq pipeline tester for task bwa
# Author: Jin Lee (leepc12@gmail.com)
import "../../chipseq.wdl" as chipseq

workflow test_bwa {
	Array[String] pe_fastqs
	Array[String] se_fastqs

	# we don't compare BAM because BAM's header includes date
	# hence md5sums don't match all the time
	String ref_pe_flagstat
	String ref_se_flagstat

	String pe_bwa_idx_tar
	String se_bwa_idx_tar

	call chipseq.bwa as pe_bwa { input :
		idx_tar = pe_bwa_idx_tar,
		fastqs = pe_fastqs,
		paired_end = true,
		cpu = 1,
	}
	call chipseq.bwa as se_bwa { input :
		idx_tar = se_bwa_idx_tar,
		fastqs = se_fastqs,
		paired_end = false,
		cpu = 1,
	}

	call chipseq.compare_md5sum { input :
		labels = [
			'pe_bwa',
			'se_bwa',
		],
		files = [
			pe_bwa.flagstat_qc,
			se_bwa.flagstat_qc,
		],
		ref_files = [
			ref_pe_flagstat,
			ref_se_flagstat,
		],
	}
}
