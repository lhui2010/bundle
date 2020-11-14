#!/usr/bin/env python

import argparse
import textwrap
import subprocess
import os
import logging

#class CTL():

#Give arguments, will return the specific CTLs
#Support Direct Print
#Support tag value change


#0 fq.gz
#1 PREFIX
#2 threads
#3 kmer_size
genomescope_sh="""
zcat {0} >{1}.cat.fq
/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/kmer/bin/jellyfish-linux count -C -m {3} -s 1000000000 -t {2} {1}.cat.fq -o {1}.jf 
/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/kmer/bin/jellyfish-linux histo -t {2} {1}.jf > {1}.reads.histo
Rscript /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/kmer/bin/genomescope/genomescope.R {1}.reads.histo {3} 150
"""

"""
bsub512 "python $BD/wrapper/kmer_wrapper.py Eu_1.fq.gz Eu_2.fq.gz "
"""

def cmd_log_and_execute(cmd):
    logging.warning(cmd)
    subprocess.run(cmd, shell = True)

def genomescope(fastq, prefix='', threads=64, kmer = 21, output= ''):
#assembly, subreads
    if(prefix == ''):
        prefix = os.path.splitext(os.path.basename(fastq[0]))[0]
    if(output == ''):
        output ="workdir_genomescope" + prefix 
    fastq_text = ' '.join(fastq)
    cmd = genomescope_sh.format(fastq_text, prefix, threads, kmer)
    cmd_log_and_execute(cmd)
    #subprocess.run(cmd, shell = True)


def main():
    prog_name = "kmer_wrapper"
    usage = "run kmer on selected fastqs"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("fastq", nargs='+',help="fastq to be evalutated in fastq.gz format")
    parser.add_argument("-t", "--threads", default=64, type=int, help="threads to run")
    args = parser.parse_args()  

    genomescope(args.fastq)
#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
