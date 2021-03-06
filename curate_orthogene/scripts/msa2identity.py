#!/usr/bin/env python
import re
import argparse
from Bio import SeqIO
import sys

#Args
parser = argparse.ArgumentParser(description='MSA To Identity ')
parser.add_argument('REF_ID', type=str, nargs = 1,
                    help='fasta of query sequences')
parser.add_argument('QRY', type=str, nargs = 1,
                    help='fasta of query sequences')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

#Prepare fasta
sys.stderr.write("Reading fasta sequences..\n")
REF_TAG = args.REF_ID[0]
qry_dict = SeqIO.to_dict(SeqIO.parse(args.QRY[0], "fasta"))
ref_dict_reverse = {}
print_buf = ""

#Loop
#First we need to determin the ref sequence used to calculate identity
#REF_TAG = "B73"
REF_LEN = 0
QRY_LEN = 0
match = {}
for qry_id in qry_dict.keys():
    if(qry_id.startswith(REF_TAG) and "-R" not in qry_id):
        REF_ID = qry_id
#        print(len(qry_dict[qry_id].seq))
#        print(len(qry_dict[qry_id].seq.__str__().replace("-", "")))
#        print(qry_dict[qry_id].seq.__str__().replace("-", ""))
        REF_LEN = len(qry_dict[qry_id].seq.__str__().replace("-", ""))
    else:
        QRY_LEN = len(qry_dict[qry_id].seq.__str__().replace("-", ""))

for qry_id in qry_dict.keys():
    match[qry_id] = 0
    if(qry_id == REF_ID):
        continue
    result_line = ""
    for pos in range(0, len(qry_dict[REF_ID].seq)):
        if(qry_dict[REF_ID].seq[pos] == "-"):
            continue
        elif(qry_dict[REF_ID].seq[pos] == qry_dict[qry_id].seq[pos]):
            match[qry_id] += 1
#        sys.stderr.write("Searching in " + qry_id + "..\n")

for k in match:
    div1 = match[k] / REF_LEN
#    div2 = match[k] / QRY_LEN
    #TODO print("{}\t{}\t{}\t{}\t{}\t{}".format(k, match[k], REF_LEN, div1, QRY_LEN, div2))
    print("{}\t{}\t{}\t{}\t{}\t{}".format(k, match[k], REF_LEN, div1, 'null', 'null'))

