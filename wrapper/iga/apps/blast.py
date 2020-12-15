"""
blast related wrappers
"""
from iga.apps.base import emain, bsub, wait_until_finish

#0 ref
#1 qry
#2 output
blastp_sh=r"""
makeblastdb -in {0} -dbtype prot
blastp -db {0} -query {1} -out {2}.bln -evalue 1e-5 -outfmt 7 -num_threads 5
"""

def blastp(ref=None, qry=None):
    r"""
    :param query:
    :param ref:
    :return:
    """
    cmd = blastp_sh.format(ref, qry, qry)
    job = bsub(cmd)
    wait_until_finish(job)
    return 0

if __name__ == "__main__":
    emain()
