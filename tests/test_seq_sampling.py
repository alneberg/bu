#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
from operator import attrgetter

from bu.seq_sampling import collect, shuffle, shatter, _shatter_recursive, sample_all
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

file_location = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(file_location,'test_data','TestE4.fa')

with open(test_data,'r') as seq_file:
    SEQS = []
    for SEQ in SeqIO.parse(seq_file, "fasta"):
        SEQS.append(SEQ)

class TestSamplerIntegration(object):
    def setup(self):
        self.seqs_l1 = sample_all(SEQS)
        
    def test_total_length(self):
        assert_true(sum(map(len, self.seqs_l1)) < sum(map(len, SEQS)),
                    msg="sampled sequences are not shorter than the originals")
    
    def test_min_length(self):
        assert_true(sum(map(len, self.seqs_l1)) > sum(map(len, SEQS))*0.5,
                    msg="sampled sequences are not shorter than the originals")

    
class TestSamplerUnit(object):
    def setup(self):
        self.seq1 = SeqRecord(Seq("A"), id = "First sequence")
        self.seq1.start_position = 0
        self.seq2 = SeqRecord(Seq("AATT"), id = "Second sequence")
        self.seq2.start_position = 0
        self.seq3 = SeqRecord(Seq("CCGG"))
        self.seq3.start_position = 0
        self.real_seq = SEQS[0]
        self.real_seq.start_position = 0

    def test_shatter(self):
        seqs1 = shatter(self.seq1, 1)
        seqs2 = []
        for i in xrange(100):
            seqs2.append(shatter(self.seq2, 1))

        
        # Keeps id
        for seq in seqs1:
            assert_equal(seq.id, self.seq1.id,
                         msg="Id not preserved")

        for seqs in seqs2:
            for seq in seqs:
                assert_equal(seq.id, self.seq2.id,
                             msg="Id not preserved")

        # Recreates the correct sequence
        for seqs in seqs2:
            seqs.sort(key=attrgetter('start_position'))
            rec_seq = ""
            for seq in seqs:
                rec_seq += str(seq.seq)
            assert_equal(rec_seq, str(self.seq2.seq),
                         msg="Sequence was not reconstructed correctly")

    def test_shuffle(self):
        seqs_l = []
        for i in xrange(100):
            seqs_l.append(shuffle([self.seq1, self.seq2, self.seq3]))

        shuffled = False
        for seqs in seqs_l:
            assert_equal(len(seqs), 3,
                         msg="Shuffled list is not of same cardinality")
            assert_true(self.seq1 in seqs,
                        msg="shuffled list is missing first sequence")
            assert_true(self.seq2 in seqs,
                        msg="shuffled list is missing second sequence")
            assert_true(self.seq3 in seqs,
                        msg="shuffled list is missing third sequence")
            if seqs != [self.seq1, self.seq2, self.seq3]:
                shuffled = True
        
        assert_true(shuffled,
                    msg = "None of 100 shufflings were shuffled")



    def test_collect(self):
        seqs1 = collect([self.seq1, self.seq2, self.seq3])
        seqs2 = collect([self.seq2, self.seq3, self.seq1])
        
        assert_true(len(seqs1) >= 3,
                    msg=("All sequences not left in the list when "
                         "last sequence length is more than 1-fraction."))

        assert_true(len(seqs1) == 3,
                    msg=("More sequences than one puts in!?"))

        assert_equal(seqs1, [self.seq1, self.seq2, self.seq3],
                     msg=("List is not preserved when last sequence "
                          "length is more than 1-fraction"))

        assert_equal(seqs2, [self.seq2, self.seq3],
                     msg=("Wrong sequence list after collect"))
        

    def test_shatter_recursive(self):
        seqs1 = _shatter_recursive(self.seq1, [], 1)
        seqs2 = _shatter_recursive(self.seq2, [], 4)
        seqs3 = _shatter_recursive(self.seq2, [self.seq3], 4)
        seqs4 = _shatter_recursive(self.seq2, [], 1)

        seqs5 = []
        for i in xrange(100):
            seqs5.append(_shatter_recursive(self.seq2, [], 1))
        
        # Stopping condition
        assert_equal(len(seqs1), 1, 
                     msg=("Single letter sequence splitted into "
                          "other than 1 sequence"))

        assert_equal(seqs2, [self.seq2],
                     msg=("Four letter sequence splitted into other "
                          "than it self when using min length 4"))

        assert_equal(seqs3, [self.seq3, self.seq2],
                     msg=("Four letter sequence splitted into other "
                          "than it self when using min length 4"))

        # Validating result
        assert_equal(sum(map(len, seqs4)), 4,
                     msg=("Resulting sequence does not conserve "
                          "sequence length."))

        assert_true(max(map(len, seqs5)) > 1,
                    msg=("Four letter sequence splitted into less "
                         "than 2 sequences every time out of 100."))

        # Start position
        assert_equal(seqs1[0].start_position, 0,
                     msg=("Single resulting sequence wrong start position"))

        for seqs in seqs5:
            for seq in seqs:
                assert_equal(str(self.seq2[seq.start_position:seq.start_position+len(seq)]), str(seq),
                             msg=("Sequence got wrong start position"))
        
