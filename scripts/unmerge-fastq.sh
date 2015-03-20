#!/bin/bash
#
#	script: unmerges interweaving paired-reads into two separate files 
#
# (1) Scenario: Paired-end forward and reverse reads interweaved by merge-paired-reads.sh
#
#          READ 1
#    |-------------->
# 5' |-------------------------------------------------| 3'
#      | | | | | | | | | | | | | | | | | | | | | | | |
# 3' |-------------------------------------------------| 5'
#                                       <--------------|
#                                            READ 2
#
#    See "Example 5: sortmerna on forward-reverse paired-end
#    reads (2 input files)" of the SortMeRNA User Manual (version 1.7 
#    and higher)
#
#    Use unmerge-paired-reads.sh to separate the interweaved reads from
#    a single file into two separate files, where READ 1 will be
#    on the same line as READ 2 but in two separate files
#
# command: bash unmerge-paired-reads.sh inputfile.fastq file1.fastq file2.fastq 
#
#
# date: May 13, 2013
# contact: evguenia.kopylova@lifl.fr
#


# check all files are given
if [ $# -lt 3 ]; then
 echo "usage: $0 merged-reads-file.fastq forward-reads-name.fastq reverse-reads-name.fastq"
 exit 2
elif [ $# -gt 3 ]; then
 echo "usage: $0 merged-reads-file.fastq forward-reads-name.fastq reverse-reads-name.fastq "
 exit 2
fi

# unmerge paired reads
echo "   Processing $2 .."
perl -pe 's/\n/\t/ if $. %4' "$1" | awk 'NR%2 {print}' | tr "\t" "\n" >| "$2"
echo "   Processing $3 .."
perl -pe 's/\n/\t/ if $. %4' "$1" | awk '(NR+1)%2 {print}' | tr "\t" "\n" >| "$3"

echo "   Done."


