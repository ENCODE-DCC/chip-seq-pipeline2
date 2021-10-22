version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bam_to_pbam {
    input {
        String dup_marker = 'picard'
        Int mapq_thresh = 30
        Boolean no_dup_removal = false
        File ref_fa
        File pe_bam
        File se_bam

        File chrsz

        File ref_pe_samtools_flagstat_qc
        File ref_se_samtools_flagstat_qc
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }
    String mito_chr_name = 'chrM'

    Int filter_cpu = 1
    Float filter_mem_factor = 0.0
    Int filter_time_hr = 24
    Float filter_disk_factor = 6.0

    call chip.filter as pe_filter { input :
        bam = pe_bam,
        no_dup_removal = false,
        paired_end = true,
        ref_fa = ref_fa,
        redact_nodup_bam = false,
        dup_marker = dup_marker,
        mapq_thresh = mapq_thresh,
        mito_chr_name = mito_chr_name,
        filter_chrs = [],
        chrsz = chrsz,

        cpu = filter_cpu,
        mem_factor = filter_mem_factor,
        picard_java_heap = '4G',
        time_hr = filter_time_hr,
        disk_factor = filter_disk_factor,
        runtime_environment = runtime_environment,
    }
    call chip.filter as se_filter { input :
        bam = se_bam,
        no_dup_removal = false,
        paired_end = false,
        ref_fa = ref_fa,
        redact_nodup_bam = false,
        dup_marker = dup_marker,
        mapq_thresh = mapq_thresh,
        mito_chr_name = mito_chr_name,
        filter_chrs = [],
        chrsz = chrsz,

        cpu = filter_cpu,
        mem_factor = filter_mem_factor,
        picard_java_heap = '4G',
        time_hr = filter_time_hr,
        disk_factor = filter_disk_factor,
        runtime_environment = runtime_environment,
    }
    call samtools_flagstat as pe_samtools_flagstat { input :
        bam = pe_filter.nodup_bam,
        runtime_environment = runtime_environment,
    }
    call samtools_flagstat as se_samtools_flagstat { input :
        bam = se_filter.nodup_bam,
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_bam_to_pbam',
            'se_bam_to_pbam',
        ],
        files = [
            pe_samtools_flagstat.flagstat_qc,
            se_samtools_flagstat.flagstat_qc,
        ],
        ref_files = [
            ref_pe_samtools_flagstat_qc,
            ref_se_samtools_flagstat_qc,
        ],
    }
}

task samtools_flagstat {
    input {
        File bam
    }
    command {
        samtools flagstat ~{bam} > flagstat.qc
    }
    output {
        File flagstat_qc = 'flagstat.qc'
    }
}