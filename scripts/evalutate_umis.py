#!/usr/bin/env python
"""A script that collects basic stats based on umi fastqs"""
import subprocess
import os
import argparse

def main(args):
    for fastq in args.input_fastqs:
        fastq_title = os.path.splitext(fastq)[0]
        nr_reads = subprocess.getoutput("gunzip -c {} | sed -n '2~4p' | wc -l".format(fastq))
        nr_unique_umis = subprocess.getoutput("gunzip -c {} | sed -n '2~4p' | sort | uniq -c | sort -n | wc -l".format(fastq))
        print(fastq_title, nr_reads, nr_unique_umis)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output_summary_stats')
    parser.add_argument('output_table')
    parser.add_argument('input_fastqs', nargs='*')
    args = parser.parse_args()
    main(args)
