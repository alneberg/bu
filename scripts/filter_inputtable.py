#!/usr/bin/env python
"""A script to select only a few ids from a concoct input table"""
import pandas as pd
from argparse import ArgumentParser
import os


def main(args):
    ids = list(pd.read_table(args.ids, header=None, squeeze=True).values)
    input_table = pd.read_table(args.input_table, index_col=0, chunksize = 10000)

    output_rows= {}
    if os.path.isfile(args.output_table):
        os.remove(args.output_table)
    with open(args.output_table, 'a') as f:
        first = True
        for chunk in input_table:
            if first:
                chunk[chunk.index.isin(ids)].to_csv(f, sep='\t', index_label = 'contig')
            else:
                chunk[chunk.index.isin(ids)].to_csv(f, sep='\t', header=None)
            first = False

if __name__=="__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('input_table', help="Concoct input table, with id as first column")
    parser.add_argument('ids', help="Ids to select, one id per row")
    parser.add_argument('output_table', help="Name of output file")
    args = parser.parse_args()
    main(args)
