#!/usr/bin/env python
"""A script to run Sarek on AWS"""
import argparse
import os,binascii
import subprocess

def random_hex():
    return binascii.b2a_hex(os.urandom(10))

def create_s3_bucket(bucket_name):
    # Untested, is at least lacking the tagging
    cmd = "aws s3api create-bucket --acl private --bucket my-bucket --region eu-west-1".format(bucket_name)
    subprocess.check_output(cmd.split(" "))

    # TODO aws s3api put-bucket-tagging

def main(args):
    # TODO Check that sample tsv exists

    # TODO Create s3 buckets if they are not supplied

    # TODO create reports dir

    # TODO run preprocessing

    # TODO run germLine variantcalling

    # TODO run annotate

    # TODO run MultiQC

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample_tsv', help="A sample tsv file for Sarek. "
                                           "Can be a local file or within an s3 bucket")
    parser.add_argument('--run_name', help="A unique name for this Sarek run. "
                                           "Will be used when naming the s3 buckets "
                                           "and the report dir. Default is todays date "
                                           "with appended random hex.", default=None)
    parser.add_argument('--genome', help="Reference iGenomes to use, default is GRCh38",
                        default='GRCh38')

    parser.add_argument('--genome_base', help="Genome base to use, default is iGenomes",
                        default="s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/{genome}")

    parser.add_argument('--results_s3', help='S3 bucket to use for results')
    parser.add_argument('--results_s3', help='S3 bucket to use for results')
    parser.add_argument('--AWS_batch_queue', help="The aws batch queue to use")
    parser.add_argument('--AWS_batch_queue_tiny', help="The tiny aws batch queue to use")

    args = parser.parse_args()
    main(args)


# Set parameters
WORKDIR="s3://sarek-work-benchmark/work-dream-ssd_2018_12_12"
OUTDIR="s3://sarek-result-benchmark/results_dream_ssd_20181212"
GENOME_BASE="s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38"
AWS_QUEUE="Sarek_pricing_benchmark_10"
AWS_QUEUE_TINY="Sarek_pricing_benchmark_10_tiny"
REPORT_DIR="Reports"

# A convenience paramete with parameters used for almost all sarek commands
COMMON_PARAMS="-profile awsbatch --awsqueue $AWS_QUEUE --awsqueue_tiny $AWS_QUEUE_TINY -work-dir $WORKDIR --outDir $OUTDIR --verbose -resume --localReportdir $REPORT_DIR"

# Check that 
