version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bwa {
    input {
        Array[File] pe_fastqs_R1
        Array[File] pe_fastqs_R2
        Array[File] se_fastqs

        # we don't compare BAM because BAM's header includes date
        # hence md5sums don't match all the time
        String ref_pe_flagstat
        String ref_se_flagstat

        String pe_bwa_idx_tar
        String se_bwa_idx_tar
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    Int bwa_cpu = 1
    Float bwa_mem_factor = 0.0
    Int bwa_time_hr = 48
    Float bwa_disk_factor = 8.0

    call chip.align as pe_bwa { input :
        aligner = 'bwa',
        idx_tar = pe_bwa_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        paired_end = true,
        use_bwa_mem_for_pe = false,
        bwa_mem_read_len_limit = 70,
        use_bowtie2_local_mode = false,
        crop_length = 0,
        crop_length_tol = 0,

        cpu = bwa_cpu,
        mem_factor = bwa_mem_factor,
        time_hr = bwa_time_hr,
        disk_factor = bwa_disk_factor,
        runtime_environment = runtime_environment,
    }
    call chip.align as se_bwa { input :
        aligner = 'bwa',
        idx_tar = se_bwa_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs,
        fastqs_R2 = [],
        paired_end = false,
        use_bwa_mem_for_pe = false,
        bwa_mem_read_len_limit = 70,
        use_bowtie2_local_mode = false,
        crop_length = 0,
        crop_length_tol = 0,

        cpu = bwa_cpu,
        mem_factor = bwa_mem_factor,
        time_hr = bwa_time_hr,
        disk_factor = bwa_disk_factor,
        runtime_environment = runtime_environment,
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
