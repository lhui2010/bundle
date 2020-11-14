#!/usr/bin/env python

import argparse
import textwrap
import subprocess

#class CTL():

#Give arguments, will return the specific CTLs
#Support Direct Print
#Support tag value change


gcpp_sh="""
pbmm2 align --sort  {0} {1} aln.ctg_id.bam
samtools index aln.ctg_id.bam
nproc={2}
gcpp --algorithm=arrow -x 5 -X 120 -q 0 -j $nproc \
        -r {0} aln.ctg_id.bam \
        -o output.fasta,output.fastq,output.vcf
"""

def flyre(fasta_file, threads=20):
#assembly, subreads
    conda_act = r"""
    export PS1="(base) \[\033]2;\h:\u $PWD\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
    source ~/lh/anaconda3/etc/profile.d/conda.sh
    conda activate flye
    flye --pacbio-corr {0} --out-dir {1} --genome-zie {2} --threads {3}
    """
    cmd = conda_act + gcpp_sh.format(ctg_file, bam_file, threads)
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
    parser.add_argument("FASTA", help="Pacbio fasta reads")
    parser.add_argument("BAM", help="subreads bam")
    parser.add_argument("-t", "--threads", default=64, type=int, help="threads")
    parser.add_argument("-o", "--output", default="workdir_fasta_name", type=str, help="workdir")
    parser.add_argument("-s", "--size", default="500m", type=str, help="size estimate of genome")
    parser.add_argument("-d", "--datatype", default="pacbio", type=str, help="input data type")
    args = parser.parse_args()  

    ctg_file = args.CONTIG
    bam_file = args.BAM
    gcpp(ctg_file, bam_file)
#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
