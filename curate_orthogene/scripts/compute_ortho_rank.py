#!/usr/bin/env python
import logging
import os
import sys
import argparse

#TODO
#Remove global variable

#OptParser
usage = """Calculate Identity from orthologs 
Usage: 
    {} Ortho_File diff.table  >diff.table.with_rank
""".format(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("ortho_file", help="The tab deliminated ortholog file of orthologs: Eg: A188_A188G12312    B73_Zm00001d001012")
parser.add_argument("diff_table", help="Diff table created by enrich_diff.py")

args = parser.parse_args()  

ortho_file = args.ortho_file
diff_table = args.diff_table

ortho_dict = {}
rank_by_gene = {}

with open(ortho_file) as fh:
    for line in fh:
        mylist = line.strip().split()
        ortho_dict[mylist[0]] = mylist[1]
        rank_by_gene[mylist[0]] = -6553465534
        rank_by_gene[mylist[1]] = 6553465534


line_number = 0
chr_multiply_factor = 1000000

with open(diff_table) as fh:
    for line in fh:
        mylist = line.rstrip('\n').split('\t')
        align_type = mylist[8]
        qry_gene = mylist[7] 
        ref_gene = mylist[12] 
        if(align_type == "="):
            line_number += 1
        for g in ([qry_gene, ref_gene]):
            if(g != "" and g in rank_by_gene):
#                logging.warning(g + "\t" + str(line_number))
                rank_by_gene[g] = line_number + chr_multiply_factor * int(mylist[0])

for qry_gene in ortho_dict:
    ref_gene = ortho_dict[qry_gene]
    assert rank_by_gene[qry_gene], "rank of {} with ref {} was not found!".format(qry_gene, ref_gene)
    assert rank_by_gene[ref_gene], "rank of {} with qry {} was not found!".format(qry_gene, ref_gene)
    rank_qry = int(rank_by_gene[qry_gene])
    rank_ref = int(rank_by_gene[ref_gene])
    diff = str(rank_ref - rank_qry)
    print("{}\t{}\t{}".format(qry_gene, ortho_dict[qry_gene], diff))
    
#1	1	2909434	2914480	A188_A188G39540-t1	0	+	A188_A188G39540-t1	<						
#1	1	5268826	5274393	B73_Zm00001d027230_T001	0	+	A188_A188G29359-t1	=	1	44289	49837	B73_Zm00001d027230_T001	0	+
#1	1	5275337	5280250	B73_Zm00001d027231_T002	0	-	A188_A188G29360-t1	=	1	50877	55716	B73_Zm00001d027231_T002	0	-
#1	1	5320587	5325291	A188_A188G29361-t1	0	+	A188_A188G29361-t1	|	1	92299	95134	B73_Zm00001d027232_T001	0	-
#1	1	5320760	5333120	A188_A188G29362-t1	0	-	A188_A188G29362-t1	|	1	111655	118312	B73_Zm00001d027233_T001	0	-
#1								>	1	118683	119739	B73_Zm00001d027234_T001	0	-
