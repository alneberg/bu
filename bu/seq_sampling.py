"""
bu/seq_sampling.py

Utilities for sampling non overlapping sequences of random length from a
set of original sequences. Useful for e.g. simulating contigs in a simplistic
fashion.
"""
from Bio import Seq
import random

def _shatter_recursive(seq, seqs, min_length):
    """Splits seq on a location that is generated from a uniform random 
    distribution. If the shortest of the resulting sequences is shorter 
    than min_length or if the longest of the resulting sequences is equal
    to the min_length, return seqs + seq, otherwise, continue splitting
    the longer sequence, and add the shorter one to seqs.
    """
    
    pos = random.randint(1, len(seq))
    s1, s2 = seq[:pos], seq[pos:]
    
    # Assign start position attribute
    s1.start_position = seq.start_position
    s2.start_position = seq.start_position + pos
    
    # Make sure s1 is of shortest length
    if len(s1) > len(s2):
        s1, s2 = s2, s1
        
    if len(s1) < min_length or len(s2) == min_length:
        seqs.append(seq)
        return seqs
    else:
        seqs.append(s1)
        return _shatter_recursive(s2, seqs, min_length)

def shatter(seq, min_length=100):
    """ keeps sequence id and start_position """
    return _shatter_recursive(seq, [], min_length)

def shuffle(seqs):
    random.shuffle(seqs)
    return seqs
    
def collect(seqs, fraction=0.8):
    """ Returns the first sequences of seqs that correspond to 
    at least fraction of the total length. 
    """
    len_seqs = map(len, seqs)
    total_len = float(sum(len_seqs))
        
    len_ack = 0
    for i, seq_l in enumerate(len_seqs):
        len_ack += seq_l
        if (len_ack / total_len) > fraction:
            return seqs[:(i+1)]
        
    return seqs

def sample_all(seqs):
    result_seqs = []
    for seq in seqs:
        seq.start_position = 0
        result_seqs += collect(shuffle(shatter(seq)))
    return result_seqs
    
