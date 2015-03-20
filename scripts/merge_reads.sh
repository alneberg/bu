#!/usr/bin/env bash

# Check all files are given
if [ $# -lt 3 ]; then
 echo "usage: $0 forward-reads.fastq.gz reverse-reads.fastq.gz merged-reads.fastq"
 exit 2
elif [ $# -gt 3 ]; then
 echo "usage: $0 forward-reads.fastq.gz reverse-reads.fastq.gz merged-reads.fastq"
 exit 2
fi

paste -d '\n' <(zcat $1 | perl -pe 's/\n/\t/ if $. %4') <(zcat $2 | perl -pe 's/\n/\t/ if $. %4')  | tr "\t" "\n" > $3
