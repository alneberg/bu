#!/usr/bin/env python
DESC="""Script to create umi fastqs from an unrelated umi fastq for several fastqs without umis.
    """

import os
from argparse import ArgumentParser
import sys
import dnaio

def main(args):
    with dnaio.open(args.umi_fastq) as umi_fh:
        umi_iter = iter(umi_fh)
        for input_fastq in args.read_fastqs:
            with dnaio.open(input_fastq) as input_fh:
                output_file = input_fastq.replace('.fastq.gz', '_fakeumi.fastq.gz')
                with dnaio.open(output_file, mode='w') as output_fh:
                    for input_record in input_fh:
                        umi_record = next(umi_iter)
                        umi_record.name = input_record.name
                        output_fh.write(umi_record)
                        

if __name__=="__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument('umi_fastq', help=("fastq file with UMI sequences but with non-matching headers"))
    parser.add_argument('read_fastqs', nargs='+', help="Fastq files from where the read name will be fetched")
    args = parser.parse_args()
    main(args)
