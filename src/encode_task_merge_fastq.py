#!/usr/bin/env python

# ENCODE DCC fastq merger wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, read_tsv, run_shell_cmd,
    strip_ext_fastq)

from encode_lib_genomic import (
    locate_trimmomatic)


def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog='ENCODE DCC fastq merger.',
                                     description='')
    parser.add_argument(
        'fastqs', nargs='+', type=str,
        help='TSV file path or list of FASTQs. '
             'FASTQs must be compressed with gzip (with .gz). '
             'Use TSV for multiple fastqs to be merged later. '
             'row=merge_id, col=end_id).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--crop-length', type=int, default=0,
                        help='Crop length for Trimmomatic.')
    parser.add_argument('--max-memory-mb', type=int,
                        help='Max memory for Trimmomatic.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # parse fastqs command line
    if args.fastqs[0].endswith('.gz') or args.fastqs[0].endswith('.fastq') or \
            args.fastqs[0].endswith('.fq'):  # it's fastq
        args.fastqs = [[f] for f in args.fastqs]  # make it a matrix
    else:  # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])

    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs) != 2:
            raise argparse.ArgumentTypeError(
                'Need 2 fastqs per replicate for paired end.')
        if not args.paired_end and len(fastqs) != 1:
            raise argparse.ArgumentTypeError(
                'Need 1 fastq per replicate for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def merge_fastqs(fastqs, end, out_dir):
    """make merged fastqs on $out_dir/R1, $out_dir/R2
    """
    out_dir = os.path.join(out_dir, end)
    mkdir_p(out_dir)
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastqs[0])))
    merged = '{}.merged.fastq.gz'.format(prefix)

    cmd = 'zcat -f {} | gzip -nc > {}'.format(
        ' '.join(fastqs),
        merged)
    run_shell_cmd(cmd)
    return merged


def trimmomatic_se(fastq, crop_len=0, nth=1, mem_mb=None):
    """crop SE fastq using trimmoatic SE mode
    """
    prefix = os.path.join(
        os.path.dirname(fastq)
        os.path.basename(strip_ext_fastq(fastq)))
    cropped = '{p}.crop_{len}.fastq.gz'.format(
        p=prefix, len=crop_len)
    assert(crop_len)
    java_mem_param = '-Xmx{}m'.format(int(mem_mb*0.9)) if mem_mb else ''

    cmd = ' '.join([
            'java', java_mem_param,
            '-jar', locate_trimmomatic(), 'SE',
            '-threads', nth,
            fastq,
            cropped,
            'MINLEN:{}'.format(crop_length),
            'CROP:{}'.format(crop_length)])
    run_shell_cmd(cmd)
    return cropped


def trimmomatic_pe(fastq_R1, fastq_R2, crop_len=0, nth=1, mem_mb=None):
    """crop PE fastq_R1 and fastq_R2 using trimmoatic PE mode
    """
    prefix_R1 = os.path.join(
        os.path.dirname(fastq_R1)
        os.path.basename(strip_ext_fastq(fastq_R1)))
    prefix_R2 = os.path.join(
        os.path.dirname(fastq_R2)
        os.path.basename(strip_ext_fastq(fastq_R2)))
    cropped_R1 = '{p}.crop_{len}.fastq.gz'.format(
        p=prefix_R1, len=crop_len)
    cropped_unpaired_R1 = '{p}.crop_{len}.unpaired.fastq.gz'.format(
        p=prefix_R1, len=crop_len)
    cropped_R2 = '{p}.crop_{len}.fastq.gz'.format(
        p=prefix_R2, len=crop_len)
    cropped_unpaired_R2 = '{p}.crop_{len}.unpaired.fastq.gz'.format(
        p=prefix_R2, len=crop_len)
    assert(crop_len)
    java_mem_param = '-Xmx{}m'.format(int(mem_mb*0.9)) if mem_mb else ''

    cmd = ' '.join([
            'java', java_mem_param,
            '-jar', locate_trimmomatic(), 'PE',
            '-threads', nth,
            fastq_R1,
            fastq_R2,
            cropped_R1,
            cropped_unpaired_R1,
            cropped_R2,
            cropped_unpaired_R2,
            'MINLEN:{}'.format(crop_length),
            'CROP:{}'.format(crop_length)])
    run_shell_cmd(cmd)
    rm_f([cropped_unpaired_R1, cropped_unpaired_R2])
    return cropped_R1, cropped_R2


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # update array with trimmed fastqs
    fastqs_R1 = []
    fastqs_R2 = []
    for fastqs in args.fastqs:
        fastqs_R1.append(fastqs[0])
        if args.paired_end:
            fastqs_R2.append(fastqs[1])

    log.info('Merging fastqs...')
    log.info('R1 to be merged: {}'.format(fastqs_R1))
    merged_R1 = merge_fastqs(fastqs_R1, 'R1', args.out_dir)
    if args.paired_end:
        log.info('R2 to be merged: {}'.format(fastqs_R2))
        merged_R2 = merge_fastqs(fastqs_R2, 'R2', args.out_dir)

    if args.crop_length:
        log.info('Cropping fastqs with crop_length {}...'.format(args.crop_length))
        if args.paired_end:
            trimmomatic_pe(merged_R1, merged_R2, args.crop_length, args.nth)
            rm_f([merged_R1, merged_R2])
        else:
            trimmomatic_se(merged_R1, args.crop_length, args.nth)
            rm_f(merged_R1)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
