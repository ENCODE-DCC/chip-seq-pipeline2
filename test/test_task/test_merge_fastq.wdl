# ENCODE DCC ChIP-Seq pipeline tester for task bowtie2
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_merge_fastq {
	# test merging rep1 and rep2
	Array[Array[String]] pe_fastqs
	Array[Array[String]] se_fastqs

	String ref_pe_merged_fastq_R1
	String ref_pe_merged_fastq_R2
	String ref_se_merged_fastq

	call chip.merge_fastq as pe_merge_fastq { input :
		fastqs = pe_fastqs,
		paired_end = true,
	}
	call chip.merge_fastq as se_merge_fastq { input :
		fastqs = se_fastqs,
		paired_end = false,
	}

	call chip.compare_md5sum { input :
		labels = [
			'pe_merge_fastq_R1',
			'pe_merge_fastq_R2',
			'se_merge_fastq',
		],
		files = [
			pe_merge_fastq.merged_fastqs[0],
			pe_merge_fastq.merged_fastqs[1],
			se_merge_fastq.merged_fastqs[0],
		],
		ref_files = [
			ref_pe_merged_fastq_R1,
			ref_pe_merged_fastq_R2,
			ref_se_merged_fastq,
		],
	}
}
