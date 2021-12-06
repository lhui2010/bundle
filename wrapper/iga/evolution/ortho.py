"""
Ortholog calculation related utils
"""
import os

from iga.apps.base import emain, bsub, sh, waitjob
import os.path as op
from iga.apps.blast import blastp, extract_top_n_hits
import logging
import re

# Require col2anc.pl located in bundle/
# 0 prefix
# 1 min_gene_in_block
# 2 max_gene_gap
# 3 qry
# 4 ref
from iga.apps.base import read_table

mcscanx_sh = """
# prepare {0}.gff and blast
format_mcscan_gff.pl {0}.gff3 > {0}.gff

MCScanX {0} -s {1} -m {2} -a
# tail -n +12 {0}.collinearity |sed "s/^#.*/###/; s/.*:\s\+//" > {0}.anchors
col2anc.pl {0}.collinearity > {0}.anchors
grep -v "#" {0}.anchors |awk '{{print $1"\t"$2}}' > {0}.ortho
awk '$3=="mRNA"'  {3}.gff3 |gff2bed.pl > {3}.bed
awk '$3=="mRNA"'  {4}.gff3 |gff2bed.pl > {4}.bed

# QRY=
# REF=
bsub  -R "span[hosts=1]" -q Q104C512G_X4  -o output.%J -e error.%J "python -m jcvi.graphics.dotplot {0}.anchors"
bsub  -R "span[hosts=1]" -q Q104C512G_X4  -o output.%J -e error.%J "python -m jcvi.compara.synteny depth --histogram {0}.anchors"
"""


def mcscanx(prefix1=None, prefix2=None, threads=4, min_gene_in_block=5, max_gene_gap=25, no_html="T", runKs='T',
            use_grid='T', top_num=10):
    """
    ⭐️️The mcscanx wrapper, execution eg:
        bsub python -m iga.evolution.ortho mcscanx Cercis_chinensis Cercis_chinensis
    :param prefix1: like Cercis_chinensis. Require Cercis_chinensis.pep and Cercis_chinensis.gff3 exists
    :param prefix2: Can be the same with prefix1, which will perform intra comparison only.
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
        if not op.exists(combine_blast + '.raw'):
            blastp(prefix1 + ".pep", prefix2 + ".pep", threads=threads, output=combine_blast + ".raw",
                   use_grid=use_grid)
        extract_top_n_hits(combine_blast + ".raw", top_num=top_num, output=combine_blast)
    sh("format_mcscan_gff.pl {0} > {1}".format(combine_gff3, combine_gff))
    # mcscanx_sh
    cmd = mcscanx_sh.format(combine_prefix, min_gene_in_block, max_gene_gap, prefix1, prefix2)
    # run
    sh(cmd)
    if runKs == 'T':
        if prefix2 == prefix1:
            combine_pep = prefix1 + '.pep'
            combine_cds = prefix1 + '.cds'
        else:
            combine_cds = combine_prefix + '.cds'
            combine_pep = combine_prefix + '.pep'
            sh("cat {0}.cds {1}.cds > {2}.cds".format(prefix1, prefix2, combine_prefix))
            sh("cat {0}.pep {1}.pep > {2}.pep".format(prefix1, prefix2, combine_prefix))
        kaks(combine_ortho, combine_cds, combine_pep, threads=threads, use_grid=use_grid, wait='F')
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
    kaks(ortho_file, cds_total, pep_total, threads=threads, use_grid=use_grid, wait='F')
    return 0


# 0 ortholog
# 1 cds fasta
# 2 pep fasta
# 3 number of threads
kaks_sh = """
#Presequitence 
#1. ParaAT (Modified to use NG86 model)
#2. KaKsCalculator
#Note:
#gene naming like following will fail
#>genea Orsat
#ParaAT will read as geneaOrsat by default

ORTHO={0}
CDS={1}
PEP={2}
echo {3} > proc
ParaAT.pl -h $ORTHO -n $CDS -a $PEP -p proc -o $ORTHO.ParaAT.out -f axt -k
join_kaks.pl $ORTHO.ParaAT.out/*.kaks >$ORTHO.kaks
rm -rf $ORTHO.ParaAT.out
"""


def kaks(ortho=None, cds=None, pep=None, threads=40, use_grid='T', wait='T'):
    """
    :param ortho: eg CORNE00006074-t2        CsaV3_1G039430
    :param cds: cds fasta
    :param pep: pep fasta
    :param use_grid: [T/F]
    :return:
    """
    cmd = kaks_sh.format(ortho, cds, pep, threads)
    if use_grid == 'T':
        jobid = bsub(cmd, name="kaks", cpus=threads)
        if wait == 'T':
            waitjob(jobid)
    else:
        sh(cmd)


def kaks_to_block(kaks=None, anchor=None):
    """
    kaks_to_block  xx-yy.kaks xx-yy.anchor
    Args:
        kaks: kaks file generated by kaks calculator
        anchor: anchor file in MCSCAN(python) format
    Returns:
        block file with ks values
    """
    # Input Example
    # ==> Lotus_japonicus.Cercis_chinensis.ortho.kaks <==
    # Sequence	Method	Ka	Ks	Ka/Ks	P-Value(Fisher)	Length	S-Sites	N-Sites	Fold-Sites(0:2:4)	Substitutions	S-Substitutions	N-Substitutions	Fold-S-Substitutions(0:2:4)	Fold-N-Substitutions(0:2:4)	Divergence-Time	Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)	GC(1:2:3)	ML-Score	AICc	Akaike-Weight	Model
    # Lj1g0000004.1.Lj1.0v1_Lojap-CECHI00000888-t1_Cechi	NG	1.47498	2.15355	0.684908	0.0163294	1824	422.829	1401.17	NA	1203	299.167	903.833	NANA	1.63228	1:1:1:1:1:1	0.374816(0.44469:0.361173:0.318584)	NA	NA	NA	NA
    # Lj1g0000004.1.Lj1.0v1_Lojap-CECHI00001701-t1_Cechi	NG	1.09922	2.17567	0.505233	2.16014e-06	1659	394.698	1264.3	NA	1009	279.75	729.25	NANA	1.35532	1:1:1:1:1:1	0.370411(0.455834:0.35169:0.303708)	NA	NA	NA	NA
    # ==> Lotus_japonicus.Cercis_chinensis.anchors <==
    # ### Alignment 0: score=812.0 e_value=1.7e-48 N=17 Contig00029_Lojap&chr05_Cechi minus
    # LjContig00029g0013512.1.Lj1.0v1_Lojap	CECHI00016458-t1_Cechi	0
    # LjContig00029g0019867.1.Lj1.0v1_Lojap	CECHI00016456-t1_Cechi	1e-58
    # LjContig00029g0017972.1.Lj1.0v1_Lojap	CECHI00016455-t1_Cechi	0
    # LjContig00029g0002447.1.Lj1.0v1_Lojap	CECHI00016448-t1_Cechi	3e-65
    # LjContig00029g0027697.1.Lj1.0v1_Lojap	CECHI00016442-t1_Cechi	8e-61
    # LjContig00029g0026581.3.Lj1.0v1_Lojap	CECHI00016436-t2_Cechi	0
    # LjContig00029g0005450.1.Lj1.0v1_Lojap	CECHI00016433-t1_Cechi	0
    # LjContig00029g0008644.1.Lj1.0v1_Lojap	CECHI00016432-t2_Cechi	3e-89
    # LjContig00029g0008860.1.Lj1.0v1_Lojap	CECHI00016427-t1_Cechi	3e-67
    genepair_to_ks = read_table(kaks)
    # with open(kaks) as fh:
    #     for line in fh:
    #         if line.startswith('Sequence'):
    #             continue
    #         mylist = line.split()
    #         genepair = mylist[0]
    #         ks = mylist[3]
    #         genepair_to_ks[genepair] = ks
    block_header = ''
    block_pair = []
    with open(anchor) as fh:
        for line in fh:
            if line.startswith('#'):
                # ### Alignment 0: score=812.0 e_value=1.7e-48 N=17 Contig00029_Lojap&chr05_Cechi minus
                if block_header != '':
                    last_block = __get_block_ks__(block_header, block_pair, genepair_to_ks)
                    print(last_block)
                block_header = line
                block_pair = []
            else:
                block_pair.append(line)
    last_block = __get_block_ks__(block_header, block_pair, genepair_to_ks)
    print(last_block)


def __get_block_ks__(block_header, block_pair, genepair_to_ks):
    """
    Inner functions
    Args:
        block_header: # ### Alignment 0: score=812.0 e_value=1.7e-48 N=17 Contig00029_Lojap&chr05_Cechi minus
        block_pair: LjContig00029g0019867.1.Lj1.0v1_Lojap	CECHI00016456-t1_Cechi	1e-58
    Returns:

    """
    sum_ks = 0
    count_ks = len(block_pair)
    new_block_pair = ""
    for p in block_pair:
        logging.debug(p)
        logging.debug(genepair_to_ks.keys())
        logging.debug(genepair_to_ks)
        mylist = p.rstrip().split()
        (qry, ref) = mylist[:2]
        pair = qry + '-' + ref
        try:
            pair_ks = genepair_to_ks[pair]['Ks']
        except KeyError:
            logging.error("Can't find gene pair {0} in Ks dictionary, continue without {0}".format(pair))
            count_ks -= 1
            continue
        new_block_pair += "\n{}\t{}\t{}".format(qry, ref, pair_ks)
        try:
            sum_ks += float(pair_ks)
        except ValueError:
            logging.error("NA value of Ks found in {}, continue without it".format(pair))
            count_ks -= 1
            continue
    mean_ks = sum_ks / count_ks
    new_block_header = block_header.rstrip() + "\t" + "Ks={}".format(mean_ks)
    return new_block_header + new_block_pair


def select_block_by_ks(anchor_ks = None, min_ks=0, max_ks=0.7):
    """
    20211206
    Args:
        anchor_ks: Input from kaks_to_block()
        min_ks:
        max_ks:

    Returns:

    """
    with open(anchor_ks) as fh:
        for line in fh:
            if line.startswith('#'):
                ks = line.rstrip().split()[-1].replace('Ks=', '')
                if ks > float(min_ks) and ks < float(max_ks):
                    flag = True
                    print(line, end = '')
            elif(flag):
                print(line, end = '')

# -----------------------------------

from iga.apps.blast import BlastTable


def rename_orthofinder_blast(seqid=None, blast=None):
    """
    rename the seqid in orthofinder blast to original gene name
    :param seqid: The "SequenceIDs.txt" produced by orthofinder
    :param blast: The Blast file (like Blast0_0.txt) produced by orthofinder
    :return: STDOUT of renamed blast file
    """
    rename_dict = {}
    with open(seqid) as fh:
        for line in fh:
            (name_abbr, name_raw) = line.rstrip().split(': ', 1)
            # Very slow, commenting the following line, change SequenceIDs.txt in advance if you don't like that name
            # name_raw = re.sub(r'\s.*', '', name_raw)
            rename_dict[name_abbr] = name_raw
    with open(blast) as fh:
        for line in fh:
            bln = BlastTable(line)
            try:
                bln.qry_id = rename_dict[bln.qry_id]
                bln.ref_id = rename_dict[bln.ref_id]
            except KeyError:
                logging.error("Can't find key {} or {}".format(bln.qry_id, bln.ref_id))
                exit(1)
            print(bln.get_line())


def mv_orthofinder_blast(SpeciesID=None, dir='.'):
    """
    rename Blast9_8.txt.out to Arabidopsis_thaliana.Solanum_penellii.blast
    :param SpeciesID:
    :param dir:
    :return:
    """
    mv_dict = {}
    with open(SpeciesID) as fh:
        for line in fh:
            (sp_abbr, name_raw) = line.rstrip().split(': ', 1)
            name_raw = name_raw.replace('.fasta', '')
            mv_dict[sp_abbr] = name_raw
    files = os.listdir(dir)
    for f in files:
        if f.startswith('Blast'):
            f_format = f.split('.')[0].replace('Blast', '')
            (qry, ref) = f_format.split('_')
            try:
                qry_new = mv_dict[qry]
                ref_new = mv_dict[ref]
            except KeyError:
                logging.error("Can't find key {} or {}".format(qry, ref))
                exit(1)
            old_name = f
            new_name = "{}.{}.blast".format(qry_new, ref_new)
            sh('mv {} {}'.format(old_name, new_name))


if __name__ == "__main__":
    emain()
