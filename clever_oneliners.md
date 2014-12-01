Clever One Liners
=================

Extract reads in fastq format that maps to a specific contig:

    samtools view -h -b mapping.bam contig_id | bamtools convert -format fastq > contig_id.fastq

Count number of sequences in a fasta file:

    grep -c '^>' seqs.fa
