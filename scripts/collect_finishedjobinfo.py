#!/usr/bin/env python
"""Script to collect finishedjobinfo statistics for
   nextflow pipelines ran on uppmax."""

import re
import os
import shutil
from argparse import ArgumentParser

DETAILS_PATH = "/sw/share/slurm/{clustername}/extended_uppmax_jobstats/{node}/{jobid}"
DETAILS_DEST_FILE = "{node}_{jobid}.stats"

def main(args):
    # Check clustername
    hostname = os.environ.get(['HOSTNAME'])
    if 'uppmax' not in hostname:
        raise Exception("This script only works on uppmax, sorry!")

    clustername_with_nr = hostname.replace('uppmax.uu.se', '')
    clustername = clustername_with_nr[:-1]
    # Fetch the slurmid from the nextflow log
    # TODO

    # run finishedjobinfo on all jobids to get the node it was ran on
    jobinfo = []
    with open(args.finishedjobinfo, 'r') as fh:
        for line in fh:
            items = line.split(' ')
            for item in items:
                # Who needs regexes, huh?
                if item.startswith('jobid'):
                    jobid = item.replace('jobid=', '')
                elif item.startswith('nodes='):
                    node = item.replace('nodes=', '')
            jobinfo.append((jobid, node))

    # fetch the extended run info for each individual job
    for jobid, node in jobinfo.items():
        src = DETAILS_PATH.format(clustername=clustername,
                                  node=node,
                                  jobid=jobid)

        dest = os.path.join(args.outdir, DETAILS_DEST_FILE.format(node=node,
                                                                  jobid=jobid))
        # copy the runinfo to a directory (since it is cleared out regularly)
        print("Copying, {}, {}".format(src, dest))
        shutil.copy(src, dest)

if __name__=="__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('finishedjobinfo', help="File with the finishedjobinfo rows that will be processed")
    parser.add_argument('outdir', help="Out directory where the stats files will be saved")
    args = parser.parse_args()
    main(args)
