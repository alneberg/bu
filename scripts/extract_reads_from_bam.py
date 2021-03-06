#!/usr/bin/env python
"""extract_reads_from_bam.py 

Extract all read pairs in a bam file where at 
least one of the mates are aligned to the references.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pysam import Samfile

def main(args):
    option = "r" if args.samformat else "rb"
    samfile = Samfile(args.bamfile, "rb")
    ref_ids = [samfile.gettid(r) for r in samfile.references]
    #Iterates over each read instead of each contig
    reads_to_print = []
    for aln in samfile.fetch(until_eof = True):
        if pair_is_aligned(aln, ref_ids):
            if args.read_pair == 1 and aln.is_read1:
                reads_to_print.append(aln)
            elif args.read_pair == 2 and aln.is_read2:
                reads_to_print.append(aln)
            elif args.read_pair == 0:
                reads_to_print.append(aln)
        if len(reads_to_print) >= 10000:
            # Flush the reads collected
            print_reads(reads_to_print)
            reads_to_print = []

    print_reads(reads_to_print)


def pair_is_aligned(aln, ref_ids):
    return aln.tid in ref_ids or aln.rnext in ref_ids

def print_reads(reads_to_print):
    # Fetch using the pyfaidx fasta object
    seq_iter = (SeqRecord(Seq(aln.seq).reverse_complement() if aln.is_reverse else Seq(aln.seq), id=aln.qname, description="") for aln in reads_to_print)
    SeqIO.write(seq_iter, sys.stdout, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-r", "--read_pair", type=int, choices=[0, 1, 2],
            default=0,
            help="Choose read pair to extract, 1 or 2, if 0 extract all.")
    parser.add_argument("bamfile", help="Input Bam file")
    parser.add_argument("-S", "--samformat", action="store_true", help="Input is in sam format, default is bam format.")

    args = parser.parse_args()

    main(args)
