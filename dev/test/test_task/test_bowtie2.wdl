version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bowtie2 {
    input {
        Array[File] pe_fastqs_R1
        Array[File] pe_fastqs_R2
        Array[File] se_fastqs

        # we don't compare BAM because BAM's header includes date
        # hence md5sums don't match all the time
        String ref_pe_flagstat
        String ref_se_flagstat    
        String ref_pe_local_flagstat
        String ref_se_local_flagstat

        String pe_bowtie2_idx_tar
        String se_bowtie2_idx_tar
    }

    Int bowtie2_cpu = 1
    Float bowtie2_mem_factor = 0.0
    Int bowtie2_time_hr = 48
    Float bowtie2_disk_factor = 5.0

    call chip.align as pe_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        paired_end = true,
        use_bwa_mem_for_pe = false,
        bwa_mem_read_len_limit = 70,
        use_bowtie2_local_mode = false,
        crop_length = 0,
        crop_length_tol = 0,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }
    call chip.align as se_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs,
        fastqs_R2 = [],
        paired_end = false,
        use_bwa_mem_for_pe = false,
        bwa_mem_read_len_limit = 70,
        use_bowtie2_local_mode = false,
        crop_length = 0,
        crop_length_tol = 0,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }

    call chip.align as pe_bowtie2_local { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        paired_end = true,
        use_bwa_mem_for_pe = false,
        bwa_mem_read_len_limit = 70,
        use_bowtie2_local_mode = true,
        crop_length = 0,
        crop_length_tol = 0,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }
    call chip.align as se_bowtie2_local { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs,
        fastqs_R2 = [],
        paired_end = false,
        use_bwa_mem_for_pe = false,
        bwa_mem_read_len_limit = 70,
        use_bowtie2_local_mode = true,
        crop_length = 0,
        crop_length_tol = 0,

        cpu = bowtie2_cpu,
        mem_factor = bowtie2_mem_factor,
        time_hr = bowtie2_time_hr,
        disk_factor = bowtie2_disk_factor,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_bowtie2',
            'se_bowtie2',
            'pe_local_bowtie2',
            'se_local_bowtie2',
        ],
        files = [
            pe_bowtie2.samstat_qc,
            se_bowtie2.samstat_qc,
            pe_bowtie2_local.samstat_qc,
            se_bowtie2_local.samstat_qc,
        ],
        ref_files = [
            ref_pe_flagstat,
            ref_se_flagstat,
            ref_pe_local_flagstat,
            ref_se_local_flagstat,
        ],
    }
}
