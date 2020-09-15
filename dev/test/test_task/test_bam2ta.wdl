version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bam2ta {
    input {
        Int bam2ta_subsample

        String pe_nodup_bam
        String se_nodup_bam

        String ref_pe_ta
        String ref_pe_ta_subsample
        String ref_se_ta
        String ref_se_ta_subsample
    }
    String mito_chr_name = 'chrM'

    Int bam2ta_cpu = 1
    Float bam2ta_mem_factor = 0.0
    Int bam2ta_time_hr = 6
    Float bam2ta_disk_factor = 4.0

    call chip.bam2ta as pe_bam2ta { input :
        bam = pe_nodup_bam,
        subsample = 0,
        paired_end = true,
        mito_chr_name = mito_chr_name,

        cpu = bam2ta_cpu,
        mem_factor = bam2ta_mem_factor,
        time_hr = bam2ta_time_hr,
        disk_factor = bam2ta_disk_factor,
    }
    call chip.bam2ta as pe_bam2ta_subsample { input :
        bam = pe_nodup_bam,
        subsample = bam2ta_subsample,
        paired_end = true,
        mito_chr_name = mito_chr_name,

        cpu = bam2ta_cpu,
        mem_factor = bam2ta_mem_factor,
        time_hr = bam2ta_time_hr,
        disk_factor = bam2ta_disk_factor,
    }
    call chip.bam2ta as se_bam2ta { input :
        bam = se_nodup_bam,
        subsample = 0,
        paired_end = false,
        mito_chr_name = mito_chr_name,

        cpu = bam2ta_cpu,
        mem_factor = bam2ta_mem_factor,
        time_hr = bam2ta_time_hr,
        disk_factor = bam2ta_disk_factor,
    }
    call chip.bam2ta as se_bam2ta_subsample { input :
        bam = se_nodup_bam,
        subsample = bam2ta_subsample,
        paired_end = false,
        mito_chr_name = mito_chr_name,

        cpu = bam2ta_cpu,
        mem_factor = bam2ta_mem_factor,
        time_hr = bam2ta_time_hr,
        disk_factor = bam2ta_disk_factor,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_bam2ta',
            'pe_bam2ta_subsample',
            'se_bam2ta',
            'se_bam2ta_subsample',
        ],
        files = [
            pe_bam2ta.ta,
            pe_bam2ta_subsample.ta,
            se_bam2ta.ta,
            se_bam2ta_subsample.ta,
        ],
        ref_files = [
            ref_pe_ta,
            ref_pe_ta_subsample,
            ref_se_ta,
            ref_se_ta_subsample,
        ],
    }
}
