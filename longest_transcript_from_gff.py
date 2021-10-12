#!/usr/bin/env python

#Extract longest transcript of a gene from GFF using GFF parser
#Author: Hui Liu
#lhui2010@gmail.com

#pip install bcbio-gff

from BCBio import GFF
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord
import sys
import re
import os
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

                length=0

                for sub_index2 in range(len(feat.sub_features[sub_index].sub_features)):

                    if(feat.sub_features[sub_index].sub_features[sub_index2].type == 'CDS'):

                        length += feat.sub_features[sub_index].sub_features[sub_index2].location.end.position - feat.sub_features[sub_index].sub_features[sub_index2].location.start.position
                if ( not(len_dict.__contains__(gene_name)) or length > len_dict[gene_name]):

                    len_dict[gene_name] = length

                    id_dict[gene_name] = feat.sub_features[sub_index].id
in_handle.close()

with open(args.in_file + '.longest_gene', 'w') as fh:
    for gene in list(id_dict):
        fh.write(id_dict[gene] + "\n")
        #fh.write(id_dict[gene], "\t", gene)
# os.system("grep -f {0}.longest_gene  {0} > {0}.longest".format(args.in_file))
