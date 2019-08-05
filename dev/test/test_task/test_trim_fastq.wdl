# ENCODE DCC ChIP-Seq pipeline tester for task bowtie2
# Author: Jin Lee (leepc12@gmail.com)
import "../../../chip.wdl" as chip
import "compare_md5sum.wdl" as compare_md5sum

workflow test_trim_fastq {
	Int xcor_pe_trim_bp = 50 				# for cross-correlation analysis only
	String pe_fastq_rep1_R1

	String ref_pe_trimmed_fastq_rep1_R1

	call chip.trim_fastq as pe_trim_fastq { input :
		fastq = pe_fastq_rep1_R1,
		trim_bp = xcor_pe_trim_bp,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'pe_trim_fastq',
		],
		files = [
			pe_trim_fastq.trimmed_fastq,
		],
		ref_files = [
			ref_pe_trimmed_fastq_rep1_R1,
		],
	}
}
