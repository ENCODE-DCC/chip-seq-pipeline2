version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_subsample_ctl {
    input {
        File pe_ta
        File se_ta
        Int pe_subsample
        Int se_subsample
        File ref_pe_ta_subsampled
        File ref_se_ta_subsampled
        File ref_pe_ta_subsampled_trivial
        File ref_se_ta_subsampled_trivial
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }
    Float subsample_ctl_mem_factor = 0.0
    Float subsample_ctl_disk_factor = 7.5

    call chip.subsample_ctl as pe_subsample_ctl { input :
        ta = pe_ta,
        paired_end = true,
        subsample = pe_subsample,

        mem_factor = subsample_ctl_mem_factor,
        disk_factor = subsample_ctl_disk_factor,
        runtime_environment = runtime_environment,
    }    
    call chip.subsample_ctl as se_subsample_ctl { input :
        ta = se_ta,
        paired_end = false,
        subsample = se_subsample,

        mem_factor = subsample_ctl_mem_factor,
        disk_factor = subsample_ctl_disk_factor,
        runtime_environment = runtime_environment,
    }
    # subsample > number of reads in TA
    # output will be just a shuffled TA.
    call chip.subsample_ctl as pe_subsample_ctl_trivial { input :
        ta = pe_ta,
        paired_end = true,
        subsample = 100000000,

        mem_factor = subsample_ctl_mem_factor,
        disk_factor = subsample_ctl_disk_factor,
        runtime_environment = runtime_environment,
    }    
    call chip.subsample_ctl as se_subsample_ctl_trivial { input :
        ta = se_ta,
        paired_end = false,
        subsample = 100000000,

        mem_factor = subsample_ctl_mem_factor,
        disk_factor = subsample_ctl_disk_factor,
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_subsample_ctl',
            'se_subsample_ctl',
            'pe_subsample_ctl_trivial',
            'se_subsample_ctl_trivial',
        ],
        files = [
            pe_subsample_ctl.ta_subsampled,
            se_subsample_ctl.ta_subsampled,
            pe_subsample_ctl_trivial.ta_subsampled,
            se_subsample_ctl_trivial.ta_subsampled,
        ],
        ref_files = [
            ref_pe_ta_subsampled,
            ref_se_ta_subsampled,
            ref_pe_ta_subsampled_trivial,
            ref_se_ta_subsampled_trivial,
        ],
    }
}
