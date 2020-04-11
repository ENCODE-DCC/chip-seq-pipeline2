# ENCODE DCC ChIP-Seq pipeline tester for task bwa
# Author: Jin Lee (leepc12@gmail.com)
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bowtie2 {
	Array[File] pe_fastqs_R1
	Array[File] pe_fastqs_R2
	Array[File] se_fastqs

	Int pe_crop_length
	Int se_crop_length
	Int pe_crop_length_tol = 2
	Int se_crop_length_tol = 2

	# we don't compare BAM because BAM's header includes date
	# hence md5sums don't match all the time
	String ref_pe_flagstat
	String ref_se_flagstat

	String ref_pe_cropped_flagstat
	String ref_se_cropped_flagstat

	String pe_bowtie2_idx_tar
	String se_bowtie2_idx_tar

	Int bowtie2_cpu = 1
	Int bowtie2_mem_mb = 20000
	Int bowtie2_time_hr = 48
	String bowtie2_disks = 'local-disk 100 HDD'

	call chip.align as pe_bowtie2 { input :
		aligner = 'bowtie2',
		idx_tar = pe_bowtie2_idx_tar,
		mito_chr_name = 'chrM',
		fastqs_R1 = pe_fastqs_R1,
		fastqs_R2 = pe_fastqs_R2,
		paired_end = true,
		use_bwa_mem_for_pe = false,
		crop_length = 0,
		crop_length_tol = 0,

		cpu = bowtie2_cpu,
		mem_mb = bowtie2_mem_mb,
		time_hr = bowtie2_time_hr,
		disks = bowtie2_disks,
	}
	call chip.align as se_bowtie2 { input :
		aligner = 'bowtie2',
		idx_tar = se_bowtie2_idx_tar,
		mito_chr_name = 'chrM',
		fastqs_R1 = se_fastqs,
		fastqs_R2 = [],
		paired_end = false,
		use_bwa_mem_for_pe = false,
		crop_length = 0,
		crop_length_tol = 0,

		cpu = bowtie2_cpu,
		mem_mb = bowtie2_mem_mb,
		time_hr = bowtie2_time_hr,
		disks = bowtie2_disks,
	}

	call chip.align as pe_cropped_bowtie2 { input :
		aligner = 'bowtie2',
		idx_tar = pe_bowtie2_idx_tar,
		mito_chr_name = 'chrM',
		fastqs_R1 = pe_fastqs_R1,
		fastqs_R2 = pe_fastqs_R2,
		paired_end = true,
		use_bwa_mem_for_pe = false,
		crop_length = pe_crop_length,
		crop_length_tol = pe_crop_length_tol,

		cpu = bowtie2_cpu,
		mem_mb = bowtie2_mem_mb,
		time_hr = bowtie2_time_hr,
		disks = bowtie2_disks,
	}
	call chip.align as se_cropped_bowtie2 { input :
		aligner = 'bowtie2',
		idx_tar = se_bowtie2_idx_tar,
		mito_chr_name = 'chrM',
		fastqs_R1 = se_fastqs,
		fastqs_R2 = [],
		paired_end = false,
		use_bwa_mem_for_pe = false,
		crop_length = se_crop_length,
		crop_length_tol = se_crop_length_tol,

		cpu = bowtie2_cpu,
		mem_mb = bowtie2_mem_mb,
		time_hr = bowtie2_time_hr,
		disks = bowtie2_disks,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'pe_bowtie2',
			'se_bowtie2',
			'pe_cropped_bowtie2',
			'se_cropped_bowtie2',
		],
		files = [
			pe_bowtie2.samstat_qc,
			se_bowtie2.samstat_qc,
			pe_cropped_bowtie2.samstat_qc,
			se_cropped_bowtie2.samstat_qc,
		],
		ref_files = [
			ref_pe_flagstat,
			ref_se_flagstat,
			ref_pe_cropped_flagstat,
			ref_se_cropped_flagstat,
		],
	}
}
