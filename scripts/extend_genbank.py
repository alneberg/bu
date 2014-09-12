#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from collections import defaultdict

def read_blast_output(blastoutfile): 
    sseq_ids = []
    records = []
    with open(blastoutfile) as in_handle:
        for line in in_handle:
            line_items = line.split("\t")
            qseq = line_items[0]
            sseq = line_items[1]
            evalue = line_items[2]
            pident = line_items[3]
            send = line_items[7]
            sstart = line_items[8]
            slen = line_items[10]

            records.append({'qseqid': qseq,
                            'sseqid': sseq,
                            'evalue': float(evalue),
                            'pident': float(pident),
                            'send': float(send),
                            'sstart': float(sstart),
                            'slen': float(slen)})
            sseq_ids.append(sseq.split('|')[2])
    return records, sseq_ids

class BlastFilterer(object):
    """Class just to be able to run blast_record_ok method within a 'filter' call
    """
    def __init__(self, scovs_threshold, evalue_threshold):
        self.scovs_threshold = scovs_threshold
        self.evalue_threshold = evalue_threshold

    def blast_record_ok(self, blast_record):
        """
        Checks for one blast record if it passes threshold filters
        """
        evalue_above_threshold = blast_record['evalue'] >= self.evalue_threshold

        alignment_length_in_subject = abs(blast_record['send'] - blast_record['sstart']) + 1
        percent_seq_covered = (alignment_length_in_subject / blast_record['slen']) * 100.0
        seq_cov_above_threshold =  percent_seq_covered >= self.scovs_threshold
    
        return evalue_above_threshold and seq_cov_above_threshold

def extend_record(rec, feature_blast_hit):
    """
    extends record rec with the feature feature_blast_hit if 
    not already present in record.
    """
    for feature in rec.features:
        if not 'locus_tag' in feature.qualifiers:
            continue
        if len(feature.qualifiers['locus_tag']) > 1:
            raise Exception("Multiple locus tags found for {0}".format(rec.id))
        feature_id = feature.qualifiers['locus_tag'][0]
        if feature_id in feature_blast_hit:
            if 'db_xref' in feature.qualifiers:
                feature.qualifiers['db_xref'].append(feature_blast_hit[feature_id])
            else:
                feature.qualifiers['db_xref'] = [feature_blast_hit[feature_id]]
    return rec

def translate_from_cdd(feature_blast_hits, cdd_all_file):
    with open(cdd_all_file, 'r') as cf:
         cddid_d = dict([(row.split('\t')[0], row.split('\t')[1].strip()) for row in cf.readlines()])
    return dict([(k, cddid_d[v['sseqid'].split('|')[-1]]) for k,v in feature_blast_hits.iteritems()])

def main(blastoutfile, genbankfile, scovs_threshold, evalue_threshold, cddid_all_file):
    # Read the blast output
    blast_records, sseq_ids = read_blast_output(blastoutfile)

    # Filter on evalue and scovs-threshold
    blast_filter = BlastFilterer(scovs_threshold, evalue_threshold)
    blast_records = filter(blast_filter.blast_record_ok, blast_records)
    
    # Convert blast records to dict indexed on qseqid
    feature_blast_hits =  dict([(blast_r['qseqid'], blast_r) for blast_r in blast_records])    

    # We only want the translated cdd in the resulting genbank file
    translated_blast_hits = translate_from_cdd(feature_blast_hits, cddid_all_file)

    new_genbank = []
    for rec in SeqIO.parse(genbankfile, 'genbank'):
        # Extend genbank record with blast hits    
        new_genbank.append(extend_record(rec, translated_blast_hits))
   
    SeqIO.write(new_genbank, sys.stdout, 'genbank') 
     


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('-b', '--blastoutfile', required=True,
           help=('Output of rpsblast run, assumed to be in tabular format whith '
               'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen. '
               'The contigs ids are assumed to be recoverable by removing the last underscore '
               'and the characters following it from the qseqid column.' ))
   parser.add_argument('-g', '--gbkfile', required=True,
           help=('GenBank file to be extended on stdout.'))
   parser.add_argument('-s', '--scovs-threshold', type=float, default=50.0,
           help='Threshold covered in percent, default=50.0')
   parser.add_argument('-e', '--evalue-threshold', type=float, default=0.0,
           help='Threshold evalue, default=0.0')
   parser.add_argument('--cddid_all_file', required=True,
           help = ('A table listing all entries in CDD together with its '
           'accessions. Used to find e.g. Pfam and TIGRFAM ids given the '
           'the CDD ids. '
           'Downloaded from: ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz'))
   args = parser.parse_args()

   main(args.blastoutfile, 
        args.gbkfile, 
        args.scovs_threshold, 
        args.evalue_threshold, 
        args.cddid_all_file)
