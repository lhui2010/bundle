#!/usr/bin/env python

import argparse
import textwrap
import subprocess

#class CTL():

#Give arguments, will return the specific CTLs
#Support Direct Print
#Support tag value change
from iga.apps.base import conda_act

#0 contig
#1 bam
#2 threads
gcpp_sh="""
set -euxo pipefail
pbmm2 align --sort  {0} {1} {0}.bam
samtools index {0}.bam
nproc={2}
gcpp --algorithm=arrow -x 5 -X 120 -q 0 -j $nproc \
        -r {0} {0}.bam \
        -o {0}.polish.fasta,{0}.polish.fastq,{0}.polish.vcf
rm {0}.bam
"""

def gcpp(ctg_file, bam_file, threads=20):
#assembly, subreads

    cmd = conda_act.format('falcon') + gcpp_sh.format(ctg_file, bam_file, threads)
#    cmd = conda_act 
    subprocess.run(cmd, shell = True)


def main():
    prog_name = "Polish Pacbio Assembly with Pacbio Reads"
    usage = "Polish Pacbio Assembly with Pacbio Reads"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("CONTIG", help="Raw assembly")
    parser.add_argument("BAM", help="subreads bam")
    parser.add_argument("-t", "--threads", default=100, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    ctg_file = args.CONTIG
    bam_file = args.BAM
    gcpp(ctg_file, bam_file)
#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
