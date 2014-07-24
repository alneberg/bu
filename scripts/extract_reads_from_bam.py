#!/usr/bin/env python
"""extract_reads_from_bam.py 

Extract all sequences from a fasta file that are 
present in a bamfile. To read from the fasta file
in an efficient manner it will be indexed, if it
not already is. This will produce a .fai file
for the fasta file."""

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pyfaidx import Fasta
from pysam import Samfile

def main(args):
    seqs = Fasta(args.reads_file)
    option = "r" if args.samformat else "rb"
    samfile = Samfile(args.bamfile, "rb")

    #Iterates over each read instead of each contig
    reads_to_print = []
    for aln in samfile.fetch(until_eof = True):
        reads_to_print.append(aln.qname)      
        if len(reads_to_print) >= 10000:
            # Flush the reads collected
            print_reads(reads_to_print)
            reads_to_print = []

    print_reads(reads_to_print)

def print_reads(reads_to_print):
    # Fetch using the pyfaidx fasta object
    seq_iter = (SeqRecord(Seq(str(seqs[read])), id=seqs[read].name, description="") for read in reads_to_print)
    SeqIO.write(seq_iter, sys.stdout, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reads_file", help="Input Fasta file.")
    parser.add_argument("bamfile", help="Input Bam file")
    parser.add_argument("-S", "--samformat", action="store_true", help="Input is in sam format, default is bam format.")

    args = parser.parse_args()

    main(args)
