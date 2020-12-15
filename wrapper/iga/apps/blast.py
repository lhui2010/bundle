"""
blast related wrappers
"""
from collections import defaultdict

from iga.apps.base import emain, bsub, wait_until_finish, logger

#0 ref
#1 qry
#2 output
blastp_sh=r"""
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
    wait_until_finish(job)
    return 0


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
            if(bitscore > highest[qry]):
                highest[qry] = bitscore
                qry_best[qry] = ref
                qry_line = line.rstrip()
            if(bitscore > highest[ref]):
                highest[ref] = bitscore
                ref_best[ref] = qry
    for k in qry_best:
        this_best_ref = qry_best[k]
        if(ref_best[this_best_ref] == k):
            #Reciprocal best
            print(qry_line[k])


if __name__ == "__main__":
    emain()
