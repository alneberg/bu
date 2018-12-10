#!/usr/bin/env python
"""
A script to split the sample sheet index information and reverse complement the R2 index.

Useful for automating cutadapt workflow. Writes output as json.
"""
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse
import sys

def main(args):
    samplesheet = pd.read_table(args.samplesheet, sep=';', names=["Sample_id", "Sample_name", "index_combo", "Seqs_delivered", "Seqs_ordered"], index_col=0)
    samplesheet["R1_index"] = samplesheet.index_combo.apply(lambda x: x.split('-')[0])
    samplesheet["R2_index"] = samplesheet.index_combo.apply(lambda x: x.split('-')[1])
    samplesheet["R2_rev_index"] = samplesheet.R2_index.apply(lambda x: str(Seq(x, IUPAC.unambiguous_dna).reverse_complement()))
    samplesheet[['R1_index', 'R2_rev_index']].to_json(sys.stdout)


if __name__ == "__main__":
   parser = argparse.ArgumentParser(description=__doc__)
   parser.add_argument('samplesheet',
           help=('samplesheet with semicolon separated columns sample_id, sample_name, index_combo, seqs_delivered, seqs_ordered. '
               'The index_combo is used to get the indices'))
   args = parser.parse_args()

   main(args)
