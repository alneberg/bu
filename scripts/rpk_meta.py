#!/usr/bin/env python
"""rpkm_meta.py 

A script to calculate RPK (number of reads per kilobase in reference) from a mapping file.
"""

import argparse
import sys
import pysam
import pandas as pd
import numpy as np

def main(args):
    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    rpks = np.zeros(len(samfile.references))    
    for i, (ref, ref_len) in enumerate(zip(samfile.references, samfile.lengths)):
        rpks[i] = samfile.count(reference=ref) / (1e-3*ref_len)

    s = pd.Series(rpks, samfile.references)
    s.to_csv(args.output_file, float_format="%.2f")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bamfile", help="Input BAM file.")
    parser.add_argument('-o', '--output_file',
        help=("Optional output file where sequences will be printed." 
            " Otherwise use stdout."))

    args = parser.parse_args()

    main(args)
