#!/usr/bin/env python

import argparse
import textwrap
import subprocess
import os
import logging

# class CTL():

# Give arguments, will return the specific CTLs
# Support Direct Print
# Support tag value change


# 0 fq.gz
# 1 PREFIX
# 2 threads
# 3 kmer_size
# 4 output
from iga.apps.base import emain, sh, waitjob, bsub, abspath_list

# 0 fastq.gz
# 1 prefix,
# 2 threads,
# 3 kmer,
# 4 output
genomescope_sh = """
zcat {0} >{1}.cat.fq
/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/kmer/bin/jellyfish-linux count -C -m {3} -s 1000000000 -t {2} {1}.cat.fq -o {1}.jf 
/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/kmer/bin/jellyfish-linux histo -t {2} {1}.jf > {1}.reads.histo
rm {1}.cat.fq {1}.jf
Rscript /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/kmer/bin/genomescope/genomescope.R {1}.reads.histo {3} 150 {4}
"""

"""
bsub512 "python $BD/wrapper/kmer_wrapper.py Eu_1.fq.gz Eu_2.fq.gz "
"""


def genomescope(fastq=None, prefix='', threads=64, kmer=21, output=''):
    """
    Estimate genome size, heterozygosity, repeat% with jellyfish and genomescope
    :param fastq:
    :param prefix:
    :param threads:
    :param kmer:
    :param output:
    :return:
    """
    # assembly, subreads
    if (prefix == ''):
        prefix = os.path.splitext(os.path.basename(fastq[0]))[0]
    if (output == ''):
        output = "workdir_genomescope" + prefix
    fastq_text = ' '.join(fastq)
    cmd = genomescope_sh.format(fastq_text, prefix, threads, kmer, output)
    bsub(cmd)
    # subprocess.run(cmd, shell = True)


# gzipped fq input is OK
# Require gce and kmer_freq in path
# 0 fq.gz
# 1 workdir
# 2 prefix
# 3 threads
# 4 kmer size
gce_sh = """

mkdir -p {1} && cd {1}

#export LD_LIBRARY_PATH=""

echo -e "Species\tGenomeSize\tHeterozygosity\tRepeat%"

ls {0} > {2}.fq.lst
kmer_freq_hash -t {3} -k {4} -l {2}.fq.lst -p {2} 2>{2}.kmerfreq.log
UNIQKMERNUM=`tail -n 11 {2}.kmerfreq.log  |head -1 |awk '{{print $2}}'`
DEPTH=`tail -n 11 {2}.kmerfreq.log  |head -1 |awk '{{print $5}}'`
gce -f {2}.freq.stat -g $UNIQKMERNUM -H 1 -c $DEPTH -b 1 >{2}.gce.out 2>{2}.gce.err

GenomeSize=`tail -n 2 {2}.gce.err |head -1 |awk '{{print $6}}'`
Heterozygosity=`tail -n 2 {2}.gce.err |head -1 |awk '{{print $7/(2-$7)/{4}}}'`
Repeat=`tail -n 2 {2}.gce.err |head -1 |awk '{{print 1-$9-$10}}'`
echo -e "{2}\t$GenomeSize\t$Heterozygosity\t$Repeat"

"""


def gce(fastq=None, prefix='', threads=64, kmer=23, workdir=''):
    """
    Estimate genome size, heterozygosity, repeat% with kmerfreq and gce
    :param fastq:
    :param prefix:
    :param threads:
    :param kmer:
    :param workdir:
    :return:
    """
    # assembly, subreads
    if prefix == '':
        prefix = os.path.basename(fastq[0]).split('.')[0]
    if workdir == '':
        workdir = "workdir_gce_" + prefix
    abspath_list(fastq)
    fastq_text = ' '.join(fastq)
    cmd = gce_sh.format(fastq_text, workdir, prefix, threads, kmer)
    job = bsub(cmd, cpus=threads, direct_submit=False, name=prefix)
    waitjob(job)
    return 0
    # subprocess.run(cmd, shell = T


# def main():
#     prog_name = "kmer_wrapper"
#     usage = "run kmer on selected fastqs"
#
#     parser = argparse.ArgumentParser(
#         prog=prog_name,
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description=textwrap.dedent(usage),
#         epilog="")
#     parser.add_argument("fastq", nargs='+',help="fastq to be evalutated in fastq.gz format")
#     parser.add_argument("-t", "--threads", default=64, type=int, help="threads to run")
#     args = parser.parse_args()
#
#     genomescope(args.fastq)
# #    flanking_distance = args.flanking

if __name__ == "__main__":
    emain()
