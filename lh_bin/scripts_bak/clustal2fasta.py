#!/usr/bin/env python


import sys
import os
from Bio import SeqIO
#import re

#python group2fasta.py $DIR/SingleCopyOrthogroups.txt total.pep cds

if(len(sys.argv) < 2):
    CRED = '\033[91m'
    CEND = '\033[0m'
    print(CRED + "Usage" + CEND)
    print ("\tpython " + sys.argv[0] + " culstal_format_file")
    print (CRED+"Example"+CEND)
    print ("\tpython group2fasta.py sh4.aln \n\tOutput: STDOUT")
    exit()

clustalw_file=sys.argv.pop()

#Get baseName
output_file = clustalw_file + ".fa"


#record_dict = SeqIO.to_dict(SeqIO.parse(clustalw_file, "clustal"))
with open(clustalw_file, "rU") as fh_in, open (output_file, mode="w") as fh_out:
    sequences = SeqIO.parse(fh_in, "clustal")
    count = SeqIO.write(sequences, fh_out, "fasta")

print ("Converted %i records" % count)

#        for record in SeqIO.parse(handle, "clustal"):
#            SeqIO.write(record, fh_out, "fasta")
            #print(record.id)
