#!/usr/bin/env python

# ENCODE DCC spp call peak wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common import *
from encode_blacklist_filter import blacklist_filter
from encode_frip import frip_shifted

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC spp callpeak',
                                        description='')
    parser.add_argument('ta', type=str,
                        help='Path for experiment IP TAGALIGN file.')
    parser.add_argument('--ctl-ta', type=str, required=True,
                        help='Path for control TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--cap-num-peak', default=300000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
    parser.add_argument('--blacklist', type=str, required=True,
                        help='Blacklist BED file.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
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

def spp(ta, ctl_ta, fraglen, cap_num_peak, nth, out_dir):
    basename_ta = os.path.basename(strip_ext_ta(ta))
    basename_ctl_ta = os.path.basename(strip_ext_ta(ctl_ta))
    basename_prefix = '{}_x_{}'.format(basename_ta, basename_ctl_ta)
    if len(basename_prefix) > 200: # UNIX cannot have filename > 255
        basename_prefix = '{}_x_control'.format(basename_ta)
    prefix = os.path.join(out_dir, basename_prefix)
    rpeak = '{}.{}.regionPeak.gz'.format(
        prefix,
        human_readable_number(cap_num_peak))
    rpeak_tmp = '{}.tmp'.format(rpeak)

    cmd0 = 'Rscript $(which run_spp.R) -c={} -i={} '
    cmd0 += '-npeak={} -odir={} -speak={} -savr={} -rf -p={}'
    cmd0 = cmd0.format(
        ta,
        ctl_ta,
        cap_num_peak,
        os.path.abspath(out_dir),
        fraglen,
        rpeak_tmp,
        nth)
    run_shell_cmd(cmd0)

    # if we have scientific representation of chr coord. then convert it to int
    cmd1 = 'zcat -f {} | awk \'BEGIN{{OFS="\\t"}}'
    cmd1 += '{{if ($2<0) $2=0; print $1,int($2),int($3),$4,$5,$6,$7,$8,$9,$10;}}\' | '
    cmd1 += 'gzip -f -nc > {}'
    cmd1 = cmd1.format(
        rpeak_tmp,
        rpeak)
    run_shell_cmd(cmd1)
    rm_f(rpeak_tmp)

    return rpeak

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks and generating signal tracks with spp...')
    rpeak = spp(args.ta, args.ctl_ta, 
        args.fraglen, args.cap_num_peak, args.nth, args.out_dir)

    if args.blacklist:
        log.info('Blacklist-filtering peaks...')
        bfilt_rpeak = blacklist_filter(
                rpeak, args.blacklist, False, args.out_dir)
    else:
        bfilt_rpeak = rpeak

    if args.ta: # if TAG-ALIGN is given
        log.info('Shifted FRiP with fragment length...')
        frip_qc = frip_shifted( args.ta, bfilt_rpeak,
            args.chrsz, args.fraglen, args.out_dir)
    else:
        frip_qc = '/dev/null'

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
