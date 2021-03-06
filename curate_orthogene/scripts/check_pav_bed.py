#!/usr/bin/env python
import re
import argparse
from Bio import SeqIO
import sys
import os
from multiprocessing import Pool

#Args
parser = argparse.ArgumentParser(description='Check whether PAV were real PAV by genblast against collinearity region')
parser.add_argument('GENE', type=str, nargs = 1,
                    help='fasta of query gene sequences')
parser.add_argument('GENOME', type=str, nargs = 1,
                    help='fasta of reference genome sequences')
parser.add_argument('BEDPE', type=str, nargs = 1,
                    help='BEDPE file of gene id and corresponding genome range')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

#TODO add argparse
threads = 50

#Prepare fasta
sys.stderr.write("Reading fasta sequences..\n")
gene_dict = SeqIO.to_dict(SeqIO.parse(args.GENE[0], "fasta"))
genome_dict = SeqIO.to_dict(SeqIO.parse(args.GENOME[0], "fasta"))
file_bedpe = args.BEDPE[0]
ref_dict_reverse = {}
print_buf = ""

root = os.getcwd()

#input bed format
#1       5259930 5344428 1       35397   119888  1       5268826 5274393 A188G29359-t1   0       +       5567    A188G29359-t1   B73_Zm00001d027230_T001
#1,2,3,10
genblast_bin = "/lustre/home/liuhui/bin/genblast_v139/genblast"
genblast_cmd_list = []
with open(file_bedpe) as fh:
    for line in fh:
        mylist = line.rstrip().split()
        for gene_id in mylist[0].rstrip(',').split(','):
            (refchr_id, start, end) = mylist[1:4]
#My bed is 1-based
            start = int(start) - 1
            if start < 0 :
                start = 0
            end = int(end)
            gene_seq = gene_dict[gene_id].seq
            genome_seq = genome_dict[refchr_id].seq[start:end]
#Write gene_seq and genome_seq
            qry_file = gene_id + ".fa"
            ref_file = gene_id + ".genome"
            out_file = qry_file + ".genblast"
            new_dir = os.path.join(root, args.BEDPE[0] + ".genblast", gene_id)

#change directory
            print("mkdir -p {}".format(new_dir))
            os.system("mkdir -p {}".format(new_dir))
            os.chdir(new_dir)
            
            with open(qry_file, "w") as fh_gene, \
                open(ref_file, "w") as fh_genome:
                fh_gene.write(">" + gene_id + "\n" + gene_seq.__str__())
                fh_genome.write(">" + "_".join([refchr_id, str(start+1), str(end)]) + "\n" + genome_seq.__str__())

            #qsub_cmd = "qsub -V -b y -cwd -N " + gene_id + " "
            genblast_cmd = "cd {} && {} -p genblastg -q {} -t {} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 10 -c 0.5 -s 0 -i 15 -x 20 -n 20 -v 2 -h 0 -j 3 -norepair -gff -cdna -pro -o {}".format(new_dir, genblast_bin, qry_file, ref_file, out_file)
            genblast_cmd_list.append(genblast_cmd)

            #os.system(qsub_cmd + genblast_cmd)

with Pool(threads) as p:
    p.map(os.system, genblast_cmd_list)

#        exit()


