"""
Ortholog calculation related utils
"""
from iga.apps.base import emain, bsub, sh
import os.path as op
from iga.apps.blast import blastp, extract_top_n_hits

# 0 prefix
mcscanx_sh="""
# prepare {0}.gff and blast
format_mcscan_gff.pl {0}.gff3 > {0}.gff

MCScanX {0} -s {1} -m {2} -a
tail -n +12 {0}.collinearity |sed "s/^#.*/###/; s/.*:\s\+//" > {0}.anchors
grep -v "#" {0}.anchors |awk '{print $1"\t"$2}' > {0}.ortho
QRY=ae
REF=ce
python -m jcvi.graphics.dotplot ${{QRY}}.${{REF}}.anchors
python -m jcvi.compara.synteny depth --histogram ${{QRY}}.${{REF}}.anchors
"""


def mcscanx(prefix1=None, prefix2=None, threads=1, min_gene_in_block=5, max_gene_gap=25, no_html="T", runKs='T',
            use_grid = 'T', top_num=10):
    """
    The mcscanx wrapper
    :param prefix1: like Cercis_chinensis. Require Cercis_chinensis.pep and Cercis_chinensis.gff3 exists
    :param prefix2: '' means self comparison
    :param min_gene_in_block:  -s 5
    :param max_gene_gap:-m 25
    :param no_html: [T/F]. T means do not produce html (MCScanX -a)
    :param comparison_type: intra (-b 1), inter (-b 2), both intra and inter (-b 0, default)
    :return:
    """
    if prefix1 == prefix2:
        combine_prefix = prefix1 + '.' + prefix1
        sh('cp {0}.gff3 {0}.{0}.gff3'.format(prefix1))
    else:
        sh('cat {0}.gff3 {1}.gff3 > {0}.{1}.gff3'.format(prefix1, prefix2))
        combine_prefix = prefix1 + '.' + prefix2
    # input
    combine_blast = combine_prefix + '.blast'
    combine_gff3 = combine_prefix + '.gff3'
    # output
    combine_gff = combine_prefix + '.gff'
    combine_ortho = combine_prefix + '.ortho'

    # prepare gff3 and blast
    if (not op.exists(combine_blast)):
        blastp(prefix1 + ".pep", prefix2 + ".pep", threads=threads, output=combine_blast + ".raw", use_grid=use_grid)
        extract_top_n_hits(combine_blast + ".raw", top_num=top_num, output=combine_blast)
    sh("format_mcscan_gff.pl {0} > {1}".format(combine_gff3, combine_gff))
    # mcscanx_sh
    cmd = mcscanx_sh.format(combine_prefix, min_gene_in_block, max_gene_gap)
    # run
    sh(cmd)
    if runKs == 'T':
        ks_dist(prefix1, prefix2, threads)
        if prefix2 == prefix1:
            combine_pep = prefix1 + '.pep'
            combine_cds = prefix1 + '.cds'
        else:
            sh("cat {0}.cds {1}.cds > {2}.cds".format(prefix1, prefix2, combine_prefix))
            sh("cat {0}.pep {1}.pep > {2}.pep".format(prefix1, prefix2, combine_prefix))
        kaks(combine_ortho, combine_pep, combine_pep, threads=threads, use_grid=use_grid)
    return 0


def ks_dist(prefix1=None, prefix2='', threads=40, use_grid='T', type='prot'):
    """
    :param prefix1: like corne
    :param prefix2: like cusat
    :param threads: 40 default
    :param use_grid: [T/F]
    :param type: [prot|nucl]
    :return: 0
    """
    if prefix2 == '':
        prefix2 = prefix1
    cmd = """
python -m jcvi.compara.catalog ortholog --dbtype {0} {1} {2}
cut -f1,2 {1}.{2}.anchors |grep -v "#" > {1}.{2}.ortho
"""
    if prefix1 == prefix2:
        cmd += """ln -s {1}.cds {1}.{2}.cds; ln -s {1}.pep {1}.{2}.pep"""
    else:
        cmd += """cat {1}.cds {2}.cds > {1}.{2}.cds; cat {1}.pep {2}.pep > {1}.{2}.pep"""
    cmd = cmd.format(type, prefix1, prefix2)
    sh(cmd)
    ortho_file = ".".join([prefix1, prefix2, 'ortho'])
    cds_total = ".".join([prefix1, prefix2, 'cds'])
    pep_total = ".".join([prefix1, prefix2, 'pep'])
    kaks(ortho_file, cds_total, pep_total, threads=threads, use_grid=use_grid)
    return 0


# 0 ortholog
# 1 cds fasta
# 2 pep fasta
# 3 number of threads
kaks_sh = """
#Presequitence 
#1. ParaAT 
#2. KaKsCalculator
#Note:
#gene naming like following will fail
#>genea Orsat
#ParaAT will read as geneaOrsat by default

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
