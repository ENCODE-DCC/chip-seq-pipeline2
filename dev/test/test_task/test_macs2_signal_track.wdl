version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_macs2_signal_track {
    input {
        Float pval_thresh

        Int fraglen
        # test macs2 for SE set only
        String se_ta
        String se_ctl_ta

        String ref_se_macs2_pval_bw # p-val signal
        String se_chrsz
        String se_gensz
    }

    Float macs2_mem_factor = 6.0
    Int macs2_time_hr = 24
    Float macs2_disk_factor = 40.0

    call chip.macs2_signal_track as se_macs2_signal_track { input :
        tas = [se_ta, se_ctl_ta],
        gensz = se_gensz,
        chrsz = se_chrsz,
        fraglen = fraglen,
        pval_thresh = pval_thresh,

        mem_factor = macs2_mem_factor,
        time_hr = macs2_time_hr,
        disk_factor = macs2_disk_factor,        
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'se_macs2_pval_bw',
        ],
        files = [
            se_macs2_signal_track.pval_bw,
        ],
        ref_files = [
            ref_se_macs2_pval_bw,
        ],
    }
}
