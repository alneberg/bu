#!/usr/bin/env python
"""split_bam_into_references.py

Split a bam file into separate bam files per 
reference with hits. """

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pysam import Samfile
from collections import defaultdict

def main(args):
    option = "r" if args.samformat else "rb"
    samfile = Samfile(args.bamfile, "rb")

    #Iterates over each read instead of each contig
    outputs = defaultdict(list)
    #import ipdb; ipdb.set_trace()
    for aln in samfile.fetch(until_eof = True):
        ref = samfile.getrname(aln.tid)
        outputs[ref].append(aln)

    for ref, alns in outputs.iteritems():
        print_reads(alns, ref, samfile.header)

def print_reads(reads_to_print, ref_name, header):
    output_name = "{0}_{1}.bam".format(args.output_base, ref_name)
    output_samfile = Samfile(output_name, "wb", header=header)
    for aln in reads_to_print:
        output_samfile.write(aln)
    output_samfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bamfile", help="Input Bam file")
    parser.add_argument("output_base", 
            help=("The base name of the output files, where "
            "the name of the reference will be appended with a "
            "underscore for each output file."))
    parser.add_argument("-S", "--samformat", action="store_true", help="Input is in sam format, default is bam format.")

    args = parser.parse_args()

    main(args)
