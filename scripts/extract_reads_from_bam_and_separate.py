#!/usr/bin/env python
"""extract_reads_from_bam_and_separate.py 

Extract all read pairs in a bam file where at 
least one of the mates are aligned to the references.
And sort these into three files depending on which mate
it is.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pysam import Samfile

def main(args):
    option = "r" if args.samformat else "rb"
    samfile = Samfile(args.bamfile, option)
    ref_ids = [samfile.gettid(r) for r in samfile.references]
    #Iterates over each read instead of each contig
    reads_to_print_1 = []
    reads_to_print_2 = []
    reads_to_print_u = []
    for aln in samfile.fetch(until_eof = True):
        if aln.tid in ref_ids: # This read is aligned
            if aln.rnext in ref_ids: # The mate is also aligned
                if aln.is_read1:
                    reads_to_print_1.append(aln)
                    reads_to_print_1 = flush_reads(reads_to_print_1, args.R1)
                elif aln.is_read2:
                    reads_to_print_2.append(aln)
                    reads_to_print_2 = flush_reads(reads_to_print_2, args.R2)
            else:
                reads_to_print_u.append(aln)
                reads_to_print_u = flush_reads(reads_to_print_u, args.u)

    print_reads(reads_to_print_1, args.R1)
    print_reads(reads_to_print_2, args.R2)
    print_reads(reads_to_print_u, args.u)

def flush_reads(reads_to_print, ofile):
    if len(reads_to_print) >= 10000:
        # Flush the reads collected
        print_reads(reads_to_print, ofile)
        return []
    else:
        return reads_to_print

def print_reads(reads_to_print, ofile):
    # Fetch using the pyfaidx fasta object
    seq_iter = (SeqRecord(Seq(aln.seq), id=aln.qname, description="") for aln in reads_to_print)
    SeqIO.write(seq_iter, open(ofile, "a"), "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bamfile", help="Input Bam file")
    parser.add_argument("-S", "--samformat", action="store_true", help="Input is in sam format, default is bam format.")
    parser.add_argument('-R1', help="Output file for all read1 with mates")
    parser.add_argument('-R2', help="Output file for all read2 with mates")
    parser.add_argument('-u', help="unpaired file for reads where the mate did not map")
    args = parser.parse_args()

    main(args)
