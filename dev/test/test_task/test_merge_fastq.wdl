# ENCODE DCC ChIP-Seq pipeline tester for task bowtie2
# Author: Jin Lee (leepc12@gmail.com)
import "../../../chip.wdl" as chip
import "compare_md5sum.wdl" as compare_md5sum

workflow test_merge_fastq {
	# test merging rep1 and rep2
	Array[File] pe_fastqs_R1
	Array[File] pe_fastqs_R2
	Array[File] se_fastqs

	File ref_pe_merged_fastq_R1
	File ref_pe_merged_fastq_R2
	File ref_se_merged_fastq

	call chip.merge_fastq as pe_merge_fastq { input :
		fastqs_R1 = pe_fastqs_R1,
		fastqs_R2 = pe_fastqs_R2,
		paired_end = true,
	}
	call chip.merge_fastq as se_merge_fastq { input :
		fastqs_R1 = se_fastqs,
		fastqs_R2 = [],
		paired_end = false,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'pe_merge_fastq_R1',
			'pe_merge_fastq_R2',
			'se_merge_fastq',
		],
		files = select_all([
			pe_merge_fastq.merged_fastq_R1,
			pe_merge_fastq.merged_fastq_R2,
			se_merge_fastq.merged_fastq_R1,
		]),
		ref_files = [
			ref_pe_merged_fastq_R1,
			ref_pe_merged_fastq_R2,
			ref_se_merged_fastq,
		],
	}
}
