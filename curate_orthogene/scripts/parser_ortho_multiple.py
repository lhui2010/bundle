#!/usr/bin/env python
import logging
import os
import sys
from multiprocessing import Pool
from optparse import OptionParser  

#TODO
#Add argparser

def main():
    usage = """Calculate Identity from orthologs 
Usage: 
    {} REF_TAG fasta A188 A188.pep >A188.rename.pep
    {} bed A188 A188.bed > A188.rename.bed
    """.format(__file__, __file__)
    parser = OptionParser(usage=usage)
      
    (options, args) = parser.parse_args()  

REF_ID = sys.argv.pop(1)
QRY_FA = sys.argv[2]
REF_FA = sys.argv[3]

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

workdir = "workdir"
threads = 55

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
def os_run(input_str, qry_fa = QRY_FA, ref_fa = REF_FA, workdir="workdir"):
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
    os.system("mkdir -p {}".format(workdir))
    
    file_lines = []
    with open(sys.argv[1]) as fh:
        file_lines = fh.readlines()

    #print(file_lines)

    qry_fa = sys.argv[2]
    ref_fa = sys.argv[3]

    with Pool(threads) as p:
        #p.apply_async(os_run, (file_lines, qry_fa, ref_fa,))
#        p.map(os_run_test, file_lines)
        p.map(os_run, file_lines)


