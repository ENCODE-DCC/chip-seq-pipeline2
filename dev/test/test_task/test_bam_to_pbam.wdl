version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bam_to_pbam {
    input {
        Array[File] pe_fastqs_R1
        Array[File] pe_fastqs_R2
        Array[File] se_fastqs
        File ref_fa

        File ref_pe_pbam_sam_gz
        File ref_se_pbam_sam_gz

        File pe_bowtie2_idx_tar
        File se_bowtie2_idx_tar
    }

    # use bowtie2 as aligner
    Int bowtie2_cpu = 1
    Float bowtie2_mem_factor = 0.0
    Int bowtie2_time_hr = 48
    Float bowtie2_disk_factor = 5.0

    call chip.align as pe_pbam_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        paired_end = true,
        use_bwa_mem_for_pe = false,
        crop_length = 0,
        crop_length_tol = 0,
        redact_bam = true,
        ref_fa = ref_fa,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }
    call chip.align as se_pbam_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs,
        fastqs_R2 = [],
        paired_end = false,
        use_bwa_mem_for_pe = false,
        crop_length = 0,
        crop_length_tol = 0,
        redact_bam = true,
        ref_fa = ref_fa,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }

    call bam_to_sam_gz as pe_bam_to_sam_gz { input:
        pbam = pe_pbam_bowtie2.bam,
    }

    call bam_to_sam_gz as se_bam_to_sam_gz { input:
        pbam = se_pbam_bowtie2.bam,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_pbam_sam_gz',
            'se_pbam_sam_gz',
        ],
        files = [
            pe_bam_to_sam_gz.sam_gz,
            se_bam_to_sam_gz.sam_gz,
        ],
        ref_files = [
            ref_pe_pbam_sam_gz,
            ref_se_pbam_sam_gz,
        ],
    }
}

task bam_to_sam_gz {
    input {
        File pbam
    }
    command {
        samtools view ~{pbam} | gzip -nc > sam.gz
    }
    output {
        File sam_gz = "sam.gz"
    }
}
