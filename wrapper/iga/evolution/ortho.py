"""
Ortholog calculation related utils
"""
from iga.apps.base import emain, bsub, sh

# 0 ortholog
# 1 cds fasta
# 2 pep fasta
# 3 number of threads
kaks_sh = """
#Presequitence 
#1. ParaAT
#2. KaKsCalculator
ORTHO={0}
CDS={1}
PEP={2}
echo {3} > proc
ParaAT.pl -h $ORTHO -n $CDS -a $PEP -p proc -o ParaAT.out -f axt -k
join_kaks.pl ParaAT.out/*.kaks >kaks_result
"""


def kaks(ortho=None, cds=None, pep=None, threads=40, use_grid='T'):
    """
    :param ortho: eg CORNE00006074-t2        CsaV3_1G039430
    :param cds: cds fasta
    :param pep: pep fasta
    :param use_grid: [T/F]
    :return:
    """
    cmd = kaks_sh.format(ortho, cds, pep, threads)
    if use_grid == 'T':
        bsub(cmd, name="kaks", cpus=threads)
    else:
        sh(cmd)


if __name__ == "__main__":
    emain()
