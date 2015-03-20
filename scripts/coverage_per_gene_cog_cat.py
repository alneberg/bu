#!/usr/bin/env python
DESC="""Script to summarize coverage over ORF:s, COGS and COG-categories.
    """

import os
import subprocess
import pandas as p
from argparse import ArgumentParser
import sys
import numpy as np

def main(args):
    df = p.read_table(args.lengths, sep='\t', names=["ORF_id", "Length"], index_col=0)
    df = df.merge( p.read_table(args.coverage, sep='\t', names=["ORF_id", "Coverage"], index_col=0), left_index=True, how='left', right_index=True)
    df = df.merge(p.read_table(args.cog_ids, sep='\t', names=["ORF_id", "COG_id"], index_col=0), left_index=True, how='left', right_index=True)
    cog_to_cat = dict((index, row.COG_Cat) for index, row in p.read_table(args.cog_to_cat, sep='\t', names=["COG_Cat", "COG_id"], index_col=1).iterrows())
    cat_to_top_cat = dict((index, (row.Desc, row.COG_Top_Cat)) for index, row in p.read_table(args.hierarchy, sep='\t', names=["COG_Cat", "Desc", "COG_Top_Cat"], index_col=0).iterrows())
    row_to_cat = {} # The intermediate level of COG categories
    row_to_header = {} # The broadest level of COG categories
    row_to_desc = {} # Description of intermediate COG Categories
    for index, row in df.iterrows():
        cog_id = row.COG_id
        if row.COG_id in cog_to_cat:
            cat = cog_to_cat[row.COG_id]
            row_to_cat[index] = cat
            row_to_desc[index] = cat_to_top_cat[cat][0]
            row_to_header[index] = cat_to_top_cat[cat][1]
    df["COG_Cat"] = p.Series(row_to_cat)
    df["Desc"] = p.Series(row_to_desc)
    df["COG_Top_Cat"] = p.Series(row_to_header)
    df["Reads"] = np.round((df.Coverage.astype(float) * df.Length.astype(float)) / 100).astype(int)
 
    df.to_csv(sys.stdout)


if __name__=="__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument('-i', '--cog-ids', help=("TSV file with ORF id and COG id respectively."))
    parser.add_argument('-m', '--cog-to-cat', help="TSV file with COG category and COG id respectively.")
    parser.add_argument('-H', '--hierarchy', help="TSV file with COG category id, description and category class respectively.")
    parser.add_argument('-c', '--coverage', help="TSV file with coverage for each ORF")
    parser.add_argument('-L', '--lengths', help="TSV file with lengths for each ORF")
    args = parser.parse_args()
    main(args)
