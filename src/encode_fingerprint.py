#!/usr/bin/env python

# ENCODE DCC fingerprint/JSD plot wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Fingerprint/JSD plot.',
                                        description='')
    parser.add_argument('bams', nargs='+', type=str,
                        help='List o paths for filtered experiment BAM files.')
    parser.add_argument('--ctl-bam', type=str,
                        help='Path for filtered control BAM file.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def bam2ta_se(bam, regex_grep_v_ta, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    ta = '{}.tagAlign.gz'.format(prefix)

    cmd = 'bedtools bamtobed -i {} | '
    cmd += 'awk \'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}\' | '
    if regex_grep_v_ta:
        cmd += 'grep -P -v \'{}\' | '.format(regex_grep_v_ta)
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        bam,
        ta)
    run_shell_cmd(cmd)
    return ta

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    log.info('Converting BAM to TAGALIGN...')
    if args.paired_end:
        ta = bam2ta_pe(args.bam, args.regex_grep_v_ta,
                        args.nth, args.out_dir)
    else:
        ta = bam2ta_se(args.bam, args.regex_grep_v_ta,
                        args.out_dir)

    if args.subsample:
        log.info('Subsampling TAGALIGN...')
        if args.paired_end:
            subsampled_ta = subsample_ta_pe(
                ta, args.subsample, False, False, args.out_dir)
        else:
            subsampled_ta = subsample_ta_se(
                ta, args.subsample, False, args.out_dir)
        temp_files.append(ta)
    else:
        subsampled_ta = ta

    if args.disable_tn5_shift:
        shifted_ta = subsampled_ta
    else:
        log.info("TN5-shifting TAGALIGN...")
        shifted_ta = tn5_shift_ta(subsampled_ta, args.out_dir)
        temp_files.append(subsampled_ta)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()