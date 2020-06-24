version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_spp {
    input {
        Int cap_num_peak

        Int fraglen
        Int ctl_subsample
        Boolean ctl_paired_end
        # test spp for SE set only
        String se_ta
        String se_ctl_ta

        String ref_se_spp_rpeak # raw narrow-peak
        String ref_se_spp_bfilt_rpeak # blacklist filtered narrow-peak
        String ref_se_spp_frip_qc
        String ref_se_spp_subsample_rpeak # raw narrow-peak
        String ref_se_spp_subsample_bfilt_rpeak # blacklist filtered narrow-peak
        String ref_se_spp_subsample_frip_qc

        String se_blacklist
        String se_chrsz
    }

    String regex_bfilt_peak_chr_name = 'chr[\\dXY]+'

    Int spp_cpu = 1
    Int spp_mem_mb = 16000
    Int spp_time_hr = 72
    String spp_disks = 'local-disk 100 HDD'

    call chip.call_peak as se_spp { input :
        peak_caller = 'spp',
        peak_type = 'regionPeak',
        gensz = se_chrsz,
        pval_thresh = 0.0,
        tas = [se_ta, se_ctl_ta],
        ctl_subsample = 0,
        ctl_paired_end = ctl_paired_end,
        chrsz = se_chrsz,
        fraglen = fraglen,
        cap_num_peak = cap_num_peak,
        blacklist = se_blacklist,
        regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name,

        cpu = spp_cpu,
        mem_mb = spp_mem_mb,
        time_hr = spp_time_hr,
        disks = spp_disks,
    }

    call chip.call_peak as se_spp_subsample { input :
        peak_caller = 'spp',
        peak_type = 'regionPeak',
        gensz = se_chrsz,
        pval_thresh = 0.0,
        tas = [se_ta, se_ctl_ta],
        ctl_subsample = ctl_subsample,
        ctl_paired_end = ctl_paired_end,
        chrsz = se_chrsz,
        fraglen = fraglen,
        cap_num_peak = cap_num_peak,
        blacklist = se_blacklist,
        regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name,

        cpu = spp_cpu,
        mem_mb = spp_mem_mb,
        time_hr = spp_time_hr,
        disks = spp_disks,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'se_spp_rpeak',
            'se_spp_bfilt_rpeak',
            'se_spp_frip_qc',
            'se_spp_subsample_rpeak',
            'se_spp_subsample_bfilt_rpeak',
            'se_spp_subsample_frip_qc',
        ],
        files = [
            se_spp.peak,
            se_spp.bfilt_peak,
            se_spp.frip_qc,
            se_spp_subsample.peak,
            se_spp_subsample.bfilt_peak,
            se_spp_subsample.frip_qc,
        ],
        ref_files = [
            ref_se_spp_rpeak,
            ref_se_spp_bfilt_rpeak,
            ref_se_spp_frip_qc,
            ref_se_spp_subsample_rpeak,
            ref_se_spp_subsample_bfilt_rpeak,
            ref_se_spp_subsample_frip_qc,
        ],
    }
}
