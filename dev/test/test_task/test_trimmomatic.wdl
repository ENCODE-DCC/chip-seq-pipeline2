version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_trimmomatic {
    input {
        Array[File] pe_fastqs_R1
        Array[File] pe_fastqs_R2
        Array[File] se_fastqs

        Int pe_crop_length
        Int se_crop_length
        Int pe_crop_length_tol = 2
        Int se_crop_length_tol = 2
        # we don't compare BAM because BAM's header includes date
        # hence md5sums don't match all the time
        File ref_pe_cropped_flagstat
        File ref_se_cropped_flagstat
        File ref_pe_cropped_phred33_flagstat
        File ref_se_cropped_phred64_flagstat

        File pe_bowtie2_idx_tar
        File se_bowtie2_idx_tar
    }

    Int bowtie2_cpu = 1
    Float bowtie2_mem_factor = 0.0
    Int bowtie2_time_hr = 48
    Float bowtie2_disk_factor = 5.0

    call chip.align as pe_cropped_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        paired_end = true,
        use_bwa_mem_for_pe = false,
        use_bowtie2_local_mode = false,
        crop_length = pe_crop_length,
        crop_length_tol = pe_crop_length_tol,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }
    call chip.align as se_cropped_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs,
        fastqs_R2 = [],
        paired_end = false,
        use_bwa_mem_for_pe = false,
        use_bowtie2_local_mode = false,
        crop_length = se_crop_length,
        crop_length_tol = se_crop_length_tol,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }

    call chip.align as pe_cropped_phred33_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        paired_end = true,
        use_bwa_mem_for_pe = false,
        use_bowtie2_local_mode = false,
        crop_length = pe_crop_length,
        crop_length_tol = pe_crop_length_tol,
        trimmomatic_phred_score_format = 'phred33',

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }
    call chip.align as se_cropped_phred64_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs,
        fastqs_R2 = [],
        paired_end = false,
        use_bwa_mem_for_pe = false,
        use_bowtie2_local_mode = false,
        crop_length = se_crop_length,
        crop_length_tol = se_crop_length_tol,
        trimmomatic_phred_score_format = 'phred64',

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_cropped_bowtie2',
            'se_cropped_bowtie2',
            'pe_cropped_phred33_bowtie2',
            'se_cropped_phred64_bowtie2',
        ],
        files = [
            pe_cropped_bowtie2.samstat_qc,
            se_cropped_bowtie2.samstat_qc,
            pe_cropped_phred33_bowtie2.samstat_qc,
            se_cropped_phred64_bowtie2.samstat_qc,
        ],
        ref_files = [
            ref_pe_cropped_flagstat,
            ref_se_cropped_flagstat,
            ref_pe_cropped_phred33_flagstat,
            ref_se_cropped_phred64_flagstat,
        ],
    }
}
