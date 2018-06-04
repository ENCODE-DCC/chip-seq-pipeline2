# ENCODE DCC ChIP-Seq pipeline tester for task bowtie2
# Author: Jin Lee (leepc12@gmail.com)
import "../../chip.wdl" as chip

workflow test_trim_fastq {
	String pe_fastq_rep1_R1

	String ref_pe_trimmed_fastq_rep1_R1

	call chip.trim_fastq as pe_trim_fastq { input :
		fastq = pe_fastq_rep1_R1,
	}

	call chip.compare_md5sum { input :
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
