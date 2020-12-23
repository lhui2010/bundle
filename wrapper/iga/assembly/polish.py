#!/usr/bin/env python

import argparse
import textwrap
import subprocess

from iga.apps.base import conda_act, get_prefix, bsub, emain, abspath_list, logger

# 0 contig.fa [abs path]
# 1 sgs.fq.gzs [abs path]
# 2 prefix
# 3 threads recommand 30
nextpolish_sh = r"""
export PATH=/ds3200_1/users_root/yitingshuang/lh/anaconda2/bin:$PATH

export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/NextPolish/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/NextPolish:$PATH

WORKDIR="workdir_nextpolish_"{2}
mkdir -p $WORKDIR
cd $WORKDIR
PREFIX={2}
contig={0}

touch sgs.fofn ; rm sgs.fofn
touch $PREFIX.0.fa ; rm $PREFIX.0.fa
     
ls {1} > sgs.fofn
ln -s $contig $PREFIX.0.fa

for i in `seq 0 1`
do
    j=`expr $i + 1`
    GENOME=$PREFIX.$i.fa
    OUTPUT=$PREFIX.$j.fa
    WORKDIR=./03_polish.sgs_round$i
    echo "[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = {3}
multithread_jobs = {3}
genome = $GENOME
genome_size = auto
workdir = $WORKDIR
polish_options = -p {{multithread_jobs}}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa
">polish_sgs.$i.cfg
#Run
    nextPolish polish_sgs.$i.cfg
    cat $WORKDIR/03.kmer_count/05.polish.ref.sh.work/polish_genome*/*.fasta > $OUTPUT
done
"""


def nextpolish(contig=None, fastq=None, threads=30):
    r"""
    :param contig: to be polished contigs
    :param fastq: fastqs, if multiple ,input with quotes like "a.fq b.fq"
    :return:
    """
    fastq_list = fastq.split()
    abspath_list([contig, fastq_list])
    logger.warning(contig)
    logger.warning(fastq_list)
    fastq = " ".join(fastq_list)
    prefix = get_prefix(contig)
    cmd = nextpolish_sh.format(contig, fastq, prefix, threads)
    bsub(cmd, name="nextpolish" + prefix, cpus=threads)


# 0 contig
# 1 bam
# 2 threads
gcpp_sh = """
set -euxo pipefail
pbmm2 align --sort  {0} {1} {0}.bam
samtools index {0}.bam
nproc={2}
gcpp --algorithm=arrow -x 5 -X 120 -q 0 -j $nproc \
        -r {0} {0}.bam \
        -o {0}.polish.fasta,{0}.polish.fastq,{0}.polish.vcf
rm {0}.bam
"""


def gcpp(contig=None, bam=None, threads=100):
    r"""
    Polish Pacbio Assembly with Pacbio Reads
    %s contig bam
    :param contig: Raw assembly
    :param bam: subreads bam
    :param threads: default=100"
    :return:
    """
    # assembly, subreads
    prefix = get_prefix(contig)
    cmd = conda_act.format('falcon') + gcpp_sh.format(contig, bam, threads)
    #    cmd = conda_act
    bsub(cmd, cpus=threads, name="gcpp" + prefix)


# def main():
#     prog_name = "Polish Pacbio Assembly with Pacbio Reads"
#     usage = "Polish Pacbio Assembly with Pacbio Reads"
#
#     parser = argparse.ArgumentParser(
#         prog=prog_name,
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description=textwrap.dedent(usage),
#         epilog="")
#     parser.add_argument("CONTIG", help="Raw assembly")
#     parser.add_argument("BAM", help="subreads bam")
#     parser.add_argument("-t", "--threads", default=100, type=int, help="flanking distance default (1000)")
#     args = parser.parse_args()
#
#     ctg_file = args.CONTIG
#     bam_file = args.BAM
#     gcpp(ctg_file, bam_file)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    emain()
