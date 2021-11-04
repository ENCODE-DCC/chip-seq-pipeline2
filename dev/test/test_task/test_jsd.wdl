version 1.0
import '../../../chip.wdl' as chip
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_jsd {
    input {
        Array[File] se_nodup_bams
        File se_ctl_nodup_bam
        File se_blacklist
        Array[File] ref_se_jsd_logs
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }
    # task level test data (BAM) is generated from BWA
    # so we keep using 30 here, this should be 255 for bowtie2 BAMs
    Int mapq_thresh = 30

    Int jsd_cpu = 1
    Float jsd_mem_factor = 0.0
    Int jsd_time_hr = 6
    Float jsd_disk_factor = 2.0

    call chip.jsd as se_jsd { input :
        nodup_bams = se_nodup_bams,
        ctl_bams = [se_ctl_nodup_bam], # use first control only
        blacklist = se_blacklist,
        mapq_thresh = mapq_thresh,

        cpu = jsd_cpu,
        mem_factor = jsd_mem_factor,
        time_hr = jsd_time_hr,
        disk_factor = jsd_disk_factor,
        runtime_environment = runtime_environment,
    }

    # take first 8 columns (vaule in other columns are random)
    scatter(i in range(2)){
        call take_8_cols { input :
            f = se_jsd.jsd_qcs[i],
        }
        call take_8_cols as ref_take_8_cols { input :
            f = ref_se_jsd_logs[i],
        }
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'se_jsd_rep1',
            'se_jsd_rep2',
        ],
        files = [
            take_8_cols.out[0],
            take_8_cols.out[1],
            #se_jsd.jsd_qcs[0],
            #se_jsd.jsd_qcs[1],
        ],
        ref_files = [
            ref_take_8_cols.out[0],
            ref_take_8_cols.out[1],
            #ref_se_jsd_logs[0],
            #ref_se_jsd_logs[1],
        ],
    }
}

task take_8_cols {
    input {
        File f
    }
    command {
        cut -f 1-8 ${f} > out.txt
    }
    output {
        File out = 'out.txt'
    }
}
