# ENCODE DCC ChIP-Seq pipeline tester for task bwa
# Author: Jin Lee (leepc12@gmail.com)
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bwa {
	Array[String] pe_fastqs
	Array[String] se_fastqs

	# we don't compare BAM because BAM's header includes date
	# hence md5sums don't match all the time
	String ref_pe_flagstat
	String ref_se_flagstat

	String pe_bwa_idx_tar
	String se_bwa_idx_tar

	Int bwa_cpu = 1
	Int bwa_mem_mb = 20000
	Int bwa_time_hr = 48
	String bwa_disks = 'local-disk 100 HDD'

	call chip.align as pe_bwa { input :
		aligner = 'bwa',
		idx_tar = pe_bwa_idx_tar,
		mito_chr_name = 'chrM',
		fastq_R1 = pe_fastqs[0],
		fastq_R2 = pe_fastqs[1],
		paired_end = true,
		use_bwa_mem_for_pe = false,

		cpu = bwa_cpu,
		mem_mb = bwa_mem_mb,
		time_hr = bwa_time_hr,
		disks = bwa_disks,
	}
	call chip.align as se_bwa { input :
		aligner = 'bwa',
		idx_tar = se_bwa_idx_tar,
		mito_chr_name = 'chrM',
		fastq_R1 = se_fastqs[0],
		paired_end = false,
		use_bwa_mem_for_pe = false,

		cpu = bwa_cpu,
		mem_mb = bwa_mem_mb,
		time_hr = bwa_time_hr,
		disks = bwa_disks,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'pe_bwa',
			'se_bwa',
		],
		files = [
			pe_bwa.samstat_qc,
			se_bwa.samstat_qc,
		],
		ref_files = [
			ref_pe_flagstat,
			ref_se_flagstat,
		],
	}
}
