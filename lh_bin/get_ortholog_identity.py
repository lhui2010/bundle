#!/usr/bin/env python
import logging
import os
import sys
from multiprocessing import Pool

#TODO
#Add argparser


#Required several scripts:
#1. /lustre/home/liuhui/bin/lh_bin/select_fasta.pl
#2. muscle in path dir
#3. msa2identity.py in current dir

#Input example
#python parser_ortho.py final_ortho.txt.add_info.format A188.pep B73.pep

def os_run_test(input_str, qry_fa = "A188.pep", ref_fa = "B73.pep", workdir="workdir"):
    print("\t".join([input_str, qry_fa, ref_fa]))

#def os_run(input_str, qry_fa, ref_fa, workdir="workdir"):
def os_run(input_str, qry_fa = "A188.pep", ref_fa = "B73.pep", workdir="workdir"):
    mylist = input_str.rstrip().split()
    qry_name = mylist[0]
    ref_name = mylist[1]
    output_fa = os.path.join(workdir, qry_name + "_" + ref_name + ".fa")
    output_aln = os.path.join(workdir, qry_name + "_" + ref_name + ".aln")
    output_identity = os.path.join(workdir, qry_name + "_" + ref_name + ".identity")
    logging.warning("{} {} {} > {} ".format("perl /lustre/home/liuhui/bin/lh_bin/select_fasta.pl", qry_name, qry_fa, output_fa))
    os.system("{} {} {} > {} ".format("perl /lustre/home/liuhui/bin/lh_bin/select_fasta.pl", qry_name, qry_fa, output_fa))
    os.system("{} {} {} >> {} ".format("perl /lustre/home/liuhui/bin/lh_bin/select_fasta.pl", ref_name, ref_fa, output_fa))
    os.system("{} {} > {} ".format("muscle -in ",  output_fa, output_aln))
    os.system("{} {} > {} ".format("python msa2identity.py ", output_aln, output_identity))

#qry_fa = "A188.pep"
#ref_fa = "B73.pep"
workdir = "workdir"
threads = 10

if __name__ == "__main__":
    os.system("mkdir -p {}".format(workdir))
    
    file_lines = []
    with open(sys.argv[1]) as fh:
        file_lines = fh.readlines()

    print(file_lines)

    qry_fa = sys.argv[2]
    ref_fa = sys.argv[3]

    with Pool(threads) as p:
        #p.apply_async(os_run, (file_lines, qry_fa, ref_fa,))
#        p.map(os_run_test, file_lines)
        p.map(os_run, file_lines)


