import sys
import os
from Bio import SeqIO
import re

#python group2fasta.py $DIR/SingleCopyOrthogroups.txt total.pep cds

if(len(sys.argv) < 4):
    CRED = '\033[91m'
    CEND = '\033[0m'
    print(CRED + "Usage" + CEND)
    print ("\tpython " + sys.argv[0] + " OrthoGroup_File FastaFile suffix")
    print (CRED+"Example"+CEND)
    print ("\tpython group2fasta.py SingleCopyOrthogroups.txt total.pep cds")
    exit()

suffix=sys.argv.pop()
fasta_file=sys.argv.pop()
group_file=sys.argv.pop()

record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

#print(record_dict["gi:12345678"])

with open (group_file )  as fh_group:
    for line in fh_group:
        genes=line.split()
        group_id = genes.pop(0)
        group_id=re.sub(":", "", group_id)
        os.system("mkdir " + group_id)
        with open(group_id + "/" + group_id + ".cds", mode="w") as fh_cds:
            for gene in genes:
                fh_cds.write(record_dict[gene].format("fasta"))
