#!/usr/bin/env python
DESC="""Script to collect finishedjobinfo statistics for
        several slurm jobs."""

import re
import os
import subprocess
import pandas as p
from argparse import ArgumentParser

def main(args):
    """For each file in slurm_dirs, parse out the name
    call finishedjobinfo for each jobid. """
    slurm_dir = args.slurm_dir
    r = re.compile("slurm-([0-9]+)\.out")
    d = {}
    for slurm_file in os.listdir(slurm_dir):
        jobids = r.findall(slurm_file)
        if len(jobids) != 1:
            continue
        jobid = jobids[0]

        info = get_info(jobid)
        d[info['jobid']] = info
        
    df = p.DataFrame.from_dict(d, orient="index")
    df.to_csv(args.output_file)

def get_info(jobid):
    p = subprocess.Popen(["finishedjobinfo", "-j", jobid], stdout=subprocess.PIPE)
    info = parse_info(p.stdout.read())
    return info

def parse_info(s):
    # Remove date before data
    data = s.split(" ")[2:]
    info = {}
    for datum in data:
        key = datum.split("=")[0]
        value = datum.split("=")[1]
        info[key] = value.strip()
    return info

if __name__=="__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument('slurm_dir', help=("Root path where all slurm output files are located, "
                                           "the names of these files will be used to parse out "
                                           "the jobid."))
    parser.add_argument('output_file', help="File where the output table will be printed.")
    args = parser.parse_args()
    main(args)
