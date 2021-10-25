"""
blast related wrappers
"""
from collections import defaultdict

from iga.apps.base import emain, bsub, waitjob, sh

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# 0 ref
# 1 qry
# 2 output
blastp_sh = r"""
makeblastdb -in {0} -dbtype prot
blastp  -num_threads 5 -db {0} -query {1} -out {2}.bln -evalue 1e-5 -outfmt 7
"""


def blastp(ref=None, qry=None, threads=5):
    r"""
    :param query:
    :param ref:
    :return:
    """
    cmd = blastp_sh.format(ref, qry, qry)
    job = bsub(cmd, cpus=threads)
    waitjob(job)
    return 0


# 0 ref
# 1 qry
# 2 output
blastn_sh = r"""
makeblastdb -in {0} -dbtype nucl
blastn  -num_threads 5 -db {0} -query {1} -out {2}.bln -evalue 1e-5 -outfmt 7
"""


def blastn(ref=None, qry=None, threads=5):
    r"""
    :param qry:
    :param ref:
    :return:
    """
    cmd = blastn_sh.format(ref, qry, qry)
    job = bsub(cmd, name='blast', cpus=5)
    # waitjob(job)
    return 0


def blast2bed(bln=None):
    """
    %s bln > bln.bed
    :param bln:
    :return:
    """
    cmd = r"""
awk '{print $1"\t"$7"\t"$8"\t"$2"\t"$12"\t."}' {0} > {0}.bed
""".format(bln)
    sh(cmd)


def filter_reciprocal_best(bln=None):
    qry_best = {}
    qry_line = {}
    ref_best = {}
    highest = defaultdict(int)
    with open(bln) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            mylist = line.rstrip().split()
            qry = mylist[0]
            ref = mylist[1]
            try:
                bitscore = float(mylist[-1])
            except ValueError:
                logger.error(line)
                continue
            if (bitscore > highest[qry]):
                highest[qry] = bitscore
                qry_best[qry] = ref
                qry_line[qry] = line.rstrip()
            if (bitscore > highest[ref]):
                highest[ref] = bitscore
                ref_best[ref] = qry
    for k in qry_best:
        this_best_ref = qry_best[k]
        if (ref_best[this_best_ref] == k):
            # Reciprocal best
            print(qry_line[k])


def extract_top_n_hits(bln=None, top_num=5):
    """
    A script function like blastall -v and -b:
    If you used to filter top 5 hits with blastall: blastall -v 5 -b 5
    You can run this script on blast file that has no filter before
    :param top_num:
    :param bln:
    :return:
    """
    cmd = """sort -k1,1 -k12,12gr -k11,11g -k3,3gr {0} > {0}.sorted.qry
# Then get the top 5 hits for every query:
for next in $(cut -f1 {0}.sorted.qry | sort -u); do grep -w -m {1} "$next" {0}.sorted.qry; done > {0}.sorted.qry.top{1}
sort -k2,2 -k12,12gr -k11,11g -k3,3gr {0} > {0}.sorted.ref
# Then get the top 5 hits for every query:
for next in $(cut -f1 {0}.sorted.ref | sort -u); do grep -w -m {1} "$next" {0}.sorted.ref; done > {0}.sorted.ref.top{1}
cat {0}.sorted.qry.top{1} {0}.sorted.ref.top{1} > {0}.top{1}
""".format(bln, top_num)
    sh(cmd)


def filter_bln(bln=None, eval=1e-5, bitscore=0):
    """
    filter blast based on evalue (1e-5 by default) and bitscore (no filter by default)
    output: STDOUT
    :param bln: the fmt6 blast file
    :param eval: float format
    :param bitscore: int format
    :return:
    """
    eval = float(eval)
    bitscore = int(bitscore)
    with open(bln) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            mylist = line.rstrip().split()
            qry = mylist[0]
            ref = mylist[1]
            try:
                this_bitscore = float(mylist[-1])
                this_eval = float(mylist[-2])
            except ValueError:
                logger.error(line)
                continue
            if bitscore == 0 or this_bitscore > bitscore:
                bitscore_true = True
            if this_eval < eval:
                eval_true = True
            if bitscore_true and eval_true:
                print(line, end='')


if __name__ == "__main__":
    emain()
