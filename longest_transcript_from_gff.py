#!/lustre/home/gaoxiaoyang/Program/python/bin/python

#Extract longest transcript of a gene from GFF using GFF parser
#Author: Hui Liu
#lhui2010@gmail.com

from BCBio import GFF
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord
import sys
#import pprint

import argparse

infile=""

parser = argparse.ArgumentParser(description='Extract longest transcripts from GFF')
parser.add_argument('in_file', nargs='?', 
                    help='GFF3 format file')
args = parser.parse_args()

if (len(sys.argv) < 2):
    parser.parse_args(['--help'])


#examiner = GFFExaminer()
#in_handle = open(in_file)
#pprint.pprint(examiner.available_limits(in_handle))
#in_handle.close()

len_dict = dict()
id_dict = dict()

#in_file = "test2.gff3"
#in_file = sys.argv[1]
#print in_file
#exit(0)
in_handle = open(args.in_file)
for rec in GFF.parse(in_handle):
    for feat in rec.features:
        gene_name = feat.id
#        print type(feat.sub_features)
#        print len(feat.sub_features)
#        continue
        if(feat.type != "gene"):
            continue
        if(len(feat.sub_features) > 0):
            for sub_index in range(len(feat.sub_features)):
                length = feat.sub_features[sub_index].location.end.position - feat.sub_features[sub_index].location.start.position
                if ( not(len_dict.has_key(gene_name)) or length > len_dict[gene_name]):
                    len_dict[gene_name] = length
                    id_dict[gene_name] = feat.sub_features[sub_index].id
            
for gene in list(id_dict):
    print id_dict[gene], "\t", gene
