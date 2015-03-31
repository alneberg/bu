#!/usr/bin/env python

DESC="""Script to handle reads mapping to OTU:s.
    """
import os
import sys
import pysam
from argparse import ArgumentParser
from collections import defaultdict, Counter

def readname(x):
    if x.is_read1:
        return "{0}_1".format(x.query_name)
    else:
        return "{0}_2".format(x.query_name)

def main(args):
    """Filter the input mapping file to only contain reads with exactly one good alignment."""
    samfile = pysam.AlignmentFile(args.mapping_file, 'rb')
    mapped_reads = defaultdict(list)

    while(1):
        try:
            r = samfile.next()
        except:
           break
        if not r.is_unmapped:
            tags = dict(r.tags)
            if 'NM' in tags:
                mapped_reads[readname(r)].append((r.reference_id, r, tags['NM']))
    
    ok_alignments = filter_alignments(mapped_reads, nm_threshold=args.nm_threshold) 
    ofile = pysam.AlignmentFile(args.output_file, 'wh', template=samfile)
    for r in ok_alignments:
        ofile.write(r)
    ofile.close()

def filter_alignments(mapped_reads, nm_threshold=3):
    """Returns alignments where exactly one alignment is better than nm_threshold.
    
    The nm_threshold is the maximum edit distance allowed from the read to the reference.
    """
    
    ok_alignments = []
    best_hit_bad = 0
    second_hit_good = 0
    for read_name, alignment_list in mapped_reads.iteritems():
        # Check if there is more than one hit for one reference
        refs = [ref_id for ref_id, r, nm in alignment_list]
        read = alignment_list[0][1]
        
        c = Counter(refs)
        max_count = c.most_common()[0][1]
        if max_count > 1:
            print "Multiple hits by one read to the same ref"
        nms = sorted([nm for ref_id, r, nm in alignment_list])
        # Make sure best hit is good enough
        if nms[0] > nm_threshold:
            best_hit_bad += 1
        
        # Make sure second best hit is not close to the best one
        elif len(nms) > 1 and nms[1] < nm_threshold:
            second_hit_good +=1
        else:
            ok_alignments.append(read)
    sys.stderr.write(("Found {0} mapping reads, where {1} had an insufficient best "
        "hit and {2} had more than one alignment that passed the threshold.\n". format(len(mapped_reads), best_hit_bad, second_hit_good)))
    return ok_alignments



if __name__=="__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument('mapping_file', help=("SAM file to read reads mapping to OTUS from."))
    parser.add_argument('output_file', help=("SAM output file."))
    parser.add_argument('--nm-threshold', default=3, help=("The maximum edit distance allowed for a mapping to be "
                        "considered good. Only reads with exactly one such alignment will be kept."))
    args = parser.parse_args()
    main(args)
