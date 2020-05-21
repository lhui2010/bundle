#!/usr/bin/env python
import logging
import os
import sys
from multiprocessing import Pool
import argparse
import subprocess

#TODO
#Remove global variable

#OptParser
usage = """Calculate Identity from orthologs 
Usage: 
    {} -o workdir REF_TAG Ortho_File QRY.fa REF.fa >A188.identity
""".format(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("REF_TAG", help="The unique keyword in gene IDs of reference genes")
parser.add_argument("OrthoFile", help="The tab deliminated ortholog file of orthologs: Eg: A188_A188G12312    B73_Zm00001d001012")
parser.add_argument("QRY_FA", help="Fasta of query fasta files")
parser.add_argument("REF_FA", help="Fasta of reference fasta files")
parser.add_argument("-o", "--output_dir", default='workdir',
                  help="specifying output directory")
parser.add_argument("-t", "--threads", default=55, type=int,
                  help="specifying threads")
args = parser.parse_args()  

WORKDIR = args.output_dir
THREADS = args.threads
REF_ID = args.REF_TAG
ORTHO_FILE = args.OrthoFile
QRY_FA = args.QRY_FA
REF_FA = args.REF_FA

#print([WORKDIR, THREADS, ORTHO_FILE, REF_ID, QRY_FA, REF_FA])
#exit()

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__)) + "/"

#Required several scripts:
#1. /lustre/home/liuhui/bin/lh_bin/select_fasta.pl
#2. muscle in path dir
#3. msa2identity.py in current dir

#Input example
#python parser_ortho.py final_ortho.txt.add_info.format A188.pep B73.pep

def os_run_test(input_str, qry_fa = "A188.pep", ref_fa = "B73.pep", workdir="workdir"):
    print("\t".join([input_str, qry_fa, ref_fa]))

#def os_run(input_str, qry_fa, ref_fa, workdir="workdir"):
#def os_run(input_str, qry_fa = "A188.pep", ref_fa = "B73.pep", workdir="workdir"):
def os_run(input_str, qry_fa = QRY_FA, ref_fa = REF_FA, workdir=WORKDIR):
    mylist = input_str.rstrip().split()
    qry_name = mylist[0]
    ref_names = mylist[1].split(',')
    output_fa = os.path.join(workdir, qry_name + ".fa")
    output_aln = os.path.join(workdir, qry_name + ".aln")
    output_identity = os.path.join(workdir, qry_name + ".identity")
    logging.warning("{} {} {} > {} ".format("perl " + SCRIPT_DIR + "select_fasta.pl", qry_name, qry_fa, output_fa))
    os.system("{} {} {} > {} ".format("perl " + SCRIPT_DIR + "select_fasta.pl", qry_name, qry_fa, output_fa))
    for ref_name in ref_names:
        logging.warning(ref_name)
        os.system("{} {} {} >> {} ".format("perl " + SCRIPT_DIR + "select_fasta.pl", ref_name, ref_fa, output_fa))
    os.system("{} {} > {} ".format("muscle -in ",  output_fa, output_aln))
    os.system("{} {} {} > {} ".format("python " + SCRIPT_DIR + "msa2identity.py ",  REF_ID, output_aln, output_identity))


if __name__ == "__main__":
    os.system("mkdir -p {}".format(WORKDIR))
    
    file_lines = []
    with open(ORTHO_FILE) as fh:
        file_lines = fh.readlines()

    #print(file_lines)

#    qry_fa = sys.argv[2]
#    ref_fa = sys.argv[3]

    with Pool(THREADS) as p:
        #p.apply_async(os_run, (file_lines, qry_fa, ref_fa,))
#        p.map(os_run_test, file_lines)
        p.map(os_run, file_lines)

    output = subprocess.check_output("for iden_ite in {}/*identity; do sort -k4,4g $iden_ite |sed -n '1p;$p' ; done".format(WORKDIR), shell=True)
    print(output.decode())
#    os.system("touch syn.identity && rm syn.identity")
#    output = os.system("for iden_ite in {}/*identity; do sort -k4,4g ${iden_ite} |sed -n '1p;$p' >>syn.identity; done".format(WORKDIR))


