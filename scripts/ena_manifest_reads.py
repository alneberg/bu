#!/usr/bin/env python
"""
Script to create manifest files for submitting raw fastqs to ENA using the webin-cli

"""

from argparse import ArgumentParser
import pandas as pd

def main(args):
    df = pd.read_table(args.input_table, sep='\t')
    for ix, row in df.iterrows():
        manifest_list = ["{} {}".format(key,value).replace('FASTQ.1', 'FASTQ') for key, value in row.to_dict().items()]
        with open("manifests/{}.txt".format(row["SAMPLE"]), 'w') as ofh:
            ofh.write("\n".join(manifest_list))  

if __name__=="__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('input_table', help=("Tab separated file with headers that will be directly translated to the manifest file"))
    args = parser.parse_args()
    main(args)
