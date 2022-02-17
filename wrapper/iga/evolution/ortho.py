"""
Ortholog calculation related utils
"""
import os

from iga.apps.base import emain, bsub, sh, waitjob, workdir_sh, mkdir
import os.path as op
from iga.apps.blast import blastp, extract_top_n_hits, extract_reciprocal_best_hits
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
            use_grid='F', top_num=10, model='YN', output=''):
    """
    Prerequisites: MCScanX(github), KaKs_Calculator, ParaAT.pl
    ⭐️️The mcscanx wrapper, execution eg:
        bsub python -m iga.evolution.ortho mcscanx Cercis_chinensis Cercis_chinensis
    :param prefix1: like Cercis_chinensis. Require Cercis_chinensis.pep and Cercis_chinensis.gff3 exists
    :param prefix2: Can be the same with prefix1, which will perform intra comparison only.
    :param min_gene_in_block:  -s 5
    :param max_gene_gap:-m 25
    :param no_html: [T/F]. T means do not produce html (MCScanX -a)
    :param comparison_type: intra (-b 1), inter (-b 2), both intra and inter (-b 0, default)
    :param model: [NG] or YN used for Ks estimation in ParaAT.pl
    :return:
    """
    combine_prefix = prefix1 + '.' + prefix2
    if output == '':
        workdir = 'work.{}'.format(combine_prefix)
    else:
        workdir = output
    #workdir_sh.format(workdir)
    mkdir(workdir)
    sh("cp {}.pep  {}/".format(prefix1, workdir))
    sh("cp {}.cds  {}/".format(prefix1, workdir))
    sh("cp {}.gff3 {}/".format(prefix1, workdir))
    sh("cp {}.pep  {}/".format(prefix2, workdir))
    sh("cp {}.cds  {}/".format(prefix2, workdir))
    sh("cp {}.gff3 {}/".format(prefix2, workdir))
    if os.path.exists(combine_prefix + ".blast"):
        sh("ln -s ../{0}.blast {1}/{0}.blast.raw".format(combine_prefix, workdir))
    os.chdir(workdir)

    if prefix1 == prefix2:
        combine_prefix = prefix1 + '.' + prefix1
        if op.exists(prefix1 + '.gff3'):
            sh('cp {0}.gff3 {0}.{0}.gff3'.format(prefix1))
        elif op.exists(prefix1 + '.gff'):
            sh('cp {0}.gff {0}.{0}.gff3'.format(prefix1))
        else:
            logging.error("You need to provide the gff3 file with format {}.gff3".format(prefix1))
            exit(1)
    else:
        if op.exists(prefix1 + '.gff3') and op.exists(prefix2 + '.gff3'):
            sh('cat {0}.gff3 {1}.gff3 > {0}.{1}.gff3'.format(prefix1, prefix2))
        elif op.exists(prefix1 + '.gff') and op.exists(prefix2 + '.gff'):
            sh('cat {0}.gff {1}.gff > {0}.{1}.gff3'.format(prefix1, prefix2))
        else:
            logging.error("You need to provide the gff3 file with format {}.gff3 {}.gff3".format(prefix1, prefix2))
            exit(1)
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
                   use_grid='F')
        extract_top_n_hits(combine_blast + ".raw", top_num=top_num, output=combine_blast, threads=threads)
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
        kaks(combine_ortho, combine_cds, combine_pep, threads=threads, use_grid=use_grid, wait='F', model=model)
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
# 4 Model
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
ParaAT.pl -h $ORTHO -n $CDS -a $PEP -p proc -o $ORTHO.ParaAT.out -f axt -k -M {4}
echo -e "Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites(0:2:4)\tSubstitutions\tS-Substitutions\tN-Substitutions\tFold-S-Substitutions(0:2:4)\tFold-N-Substitutions(0:2:4)\tDivergence-Time\tSubstitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\tGC(1:2:3)\tML-Score\tAICc\tAkaike-Weight\tModel"> $ORTHO.kaks 
find $ORTHO.ParaAT.out/ -name '*.kaks' |xargs tail -q -n 1 >>$ORTHO.kaks
# join_kaks.pl $ORTHO.ParaAT.out/*.kaks >$ORTHO.kaks
rm -rf $ORTHO.ParaAT.out
"""


def kaks(ortho=None, cds=None, pep=None, threads=40, use_grid='T', wait='T', model='YN'):
    """
    :param ortho: eg CORNE00006074-t2        CsaV3_1G039430
    :param cds: cds fasta
    :param pep: pep fasta
    :param use_grid: [T/F]
    :return:
    """
    cmd = kaks_sh.format(ortho, cds, pep, threads, model)
    if use_grid == 'T':
        jobid = bsub(cmd, name="kaks", cpus=threads, options='-m yi04')
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
        STDOUT
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
    logging.info("Reading ks table with DictReader")
    genepair_to_ks = read_table(kaks)
    logging.info("Reading ks table Complete")
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
    if last_block is not None:
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
        # logging.debug(p)
        # logging.debug(genepair_to_ks.keys())
        # logging.debug(genepair_to_ks)
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
    try:
        mean_ks = sum_ks / count_ks
    except ZeroDivisionError:
        logging.error("No ks found in block {}".format(block_header))
        return None
    new_block_header = block_header.rstrip() + "\t" + "Ks={}".format(mean_ks)
    return new_block_header + new_block_pair


def select_block_by_ks(anchor_ks=None, min_ks=0, max_ks=0.7):
    """
    20211206
    Args:
        anchor_ks: Input from kaks_to_block()
        min_ks:
        max_ks:

    Returns:

    """
    flag = False
    with open(anchor_ks) as fh:
        for line in fh:
            if line.startswith('#'):
                ks = line.rstrip().split()[-1].replace('Ks=', '')
                ks = float(ks)
                if ks > float(min_ks) and ks < float(max_ks):
                    flag = True
                    print(line, end='')
                else:
                    flag = False
            elif flag:
                print(line, end='')

# work inside a tmp directory
# 0: xx.yy.ortho
# 1: xx.xx.ortho
# 2: yy.yy.ortho
commonWGD_sh = r"""
export OMP_NUM_THREADS="8";
selectItem.pl 0 0,1 ../{0}       ../{1} |sed "s/\t/-/" > {0}.left
selectItem.pl 1 0,1 ../{0}       ../{2} |sed "s/\t/-/" > {0}.right
selectItem.pl -h    {0}.left  ../{1}.kaks >{0}.left.kaks
selectItem.pl -h    {0}.right ../{2}.kaks >{0}.right.kaks
get_ks_peak.py      {0}.left.kaks > ../{0}.left.peak
get_ks_peak.py      {0}.right.kaks > ../{0}.right.peak
get_ks_peak.py      ../{0}.kaks > ../{0}.orthopeak
"""


def commonWGD1(wgd_ortho=None):
    """
    Input Eg: Medicago_truncatula.Senna_tora.ortho
    Test if there was common WGD between xx and yy given xx.yy.ortho and xx.yy.ortho.ks as well as xx.xx.ortho[.ks] and
    yy.yy.ortho[.ks]
    Args:
        wgd_ortho: like Medicago_truncatula.Senna_tora.ortho
        requires the existence of Medicago_truncatula.Medicago_truncatula.ortho and Medicago_truncatula
        also requires Medicago_truncatula.Senna_tora.ortho.kaks
    Returns:
    """
    cmd = workdir_sh.format('work_' + wgd_ortho)
    mylist = wgd_ortho.split('.')
    wgd_left = '.'.join([mylist[0], mylist[0], 'ortho'])
    wgd_right = '.'.join([mylist[1], mylist[1], 'ortho'])
    if wgd_right == wgd_left:
        logging.info("Ignoring paralogous")
        return 0
    cmd += commonWGD_sh.format(wgd_ortho, wgd_left, wgd_right)
    bsub(cmd, cpus=8, name='GMM_peak')




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


def rnaseqks_lh(prefix1=None, prefix2=None, threads=4, runKs='T',
             use_grid='F', eval=1e-5, max_hits=10, model='YN', iden=0.2, output=''):
    r"""
    The LiuHui pipeline for calculating paralog Ks (similar to genome Ks calculation)
    Prerequisites: blast KaKs_Calculator, ParaAT.pl
    ⭐️️The rnaseqks wrapper, execution eg:
        bsub python -m iga.evolution.ortho mcscanx Cercis_chinensis Cercis_chinensis
    : param prefix1: like Cercis_chinensis. Require Cercis_chinensis.pep and Cercis_chinensis.gff3 exists
    : param prefix2: Can be the same with prefix1, which will perform intra comparison only.
    : param model: YN[default] or NG used for Ks estimation in ParaAT.pl (modified to accept model arg input)
    : param iden: min identity for blastp result (default 0.2 ), cannot be modified
    : param max_hits: 10, cannot be modified.
    :return:
    """
    combine_prefix = prefix1 + '.' + prefix2
    if output == '':
        workdir = 'work.{}'.format(combine_prefix)
    else:
        workdir = output
    #workdir_sh.format(workdir)
    mkdir(workdir)
    sh("cp {}.pep  {}/".format(prefix1, workdir))
    sh("cp {}.cds  {}/".format(prefix1, workdir))
    sh("cp {}.pep  {}/".format(prefix2, workdir))
    sh("cp {}.cds  {}/".format(prefix2, workdir))
    os.chdir(workdir)

    # input
    combine_blast = combine_prefix + '.blast'
    # output
    combine_ortho = combine_prefix + '.ortho'

    #Protein identity
    piden=20

    # prepare gff3 and blast
    if (not op.exists(combine_blast)):
        if not op.exists(combine_blast + '.raw'):
            blastp(prefix1 + ".pep", prefix2 + ".pep", threads=threads, output=combine_blast + ".raw",
                   use_grid='F', eval=eval)
        #sh("awk '$3>={2}' {0} | onlyten.pl > {1}".format(combine_blast + ".raw", combine_blast, piden))
        extract_top_n_hits(combine_blast + ".raw", top_num=max_hits, output=combine_blast, threads=threads, iden=iden)
    sh("cut -f1,2 {0} > {1}.raw".format(combine_blast, combine_ortho))
    uniq_ortho(combine_ortho + '.raw', combine_ortho)

    if runKs == 'T':
        if prefix2 == prefix1:
            combine_pep = prefix1 + '.pep'
            combine_cds = prefix1 + '.cds'
        else:
            combine_cds = combine_prefix + '.cds'
            combine_pep = combine_prefix + '.pep'
            sh("cat {0}.cds {1}.cds > {2}.cds".format(prefix1, prefix2, combine_prefix))
            sh("cat {0}.pep {1}.pep > {2}.pep".format(prefix1, prefix2, combine_prefix))
        kaks(combine_ortho, combine_cds, combine_pep, threads=threads, use_grid=use_grid, wait='F', model=model)
    return 0


def rnaseqks(prefix1=None, prefix2=None, threads=4, runKs='T',
             use_grid='F', eval=10, max_hits=10, model='YN', iden=0.2, output=''):
    r"""
    The YangYa(2015,MBE) pipeline for calculating paralog Ks
    Prerequisites: blast KaKs_Calculator, ParaAT.pl
    ⭐️️The rnaseqks wrapper, execution eg:
        bsub python -m iga.evolution.ortho mcscanx Cercis_chinensis Cercis_chinensis
    : param prefix1: like Cercis_chinensis. Require Cercis_chinensis.pep and Cercis_chinensis.gff3 exists
    : param prefix2: Can be the same with prefix1, which will perform intra comparison only.
    : param model: YN[default] or NG used for Ks estimation in ParaAT.pl (modified to accept model arg input)
    : param iden: min identity for blastp result (default 0.2 ), cannot be modified
    : param max_hits: 10, cannot be modified.
    :return:
    """
    combine_prefix = prefix1 + '.' + prefix2
    if output == '':
        workdir = 'work.{}'.format(combine_prefix)
    else:
        workdir = output
    #workdir_sh.format(workdir)
    mkdir(workdir)
    sh("cp {}.pep  {}/".format(prefix1, workdir))
    sh("cp {}.cds  {}/".format(prefix1, workdir))
    sh("cp {}.pep  {}/".format(prefix2, workdir))
    sh("cp {}.cds  {}/".format(prefix2, workdir))
    os.chdir(workdir)

    # input
    combine_blast = combine_prefix + '.blast'
    # output
    combine_ortho = combine_prefix + '.ortho'

    #Protein identity
    piden=20

    # prepare gff3 and blast
    if (not op.exists(combine_blast)):
        if not op.exists(combine_blast + '.raw'):
            blastp(prefix1 + ".pep", prefix2 + ".pep", threads=threads, output=combine_blast + ".raw",
                   use_grid='F', eval=eval)
        sh("awk '$3>={2}' {0} | onlyten.pl > {1}".format(combine_blast + ".raw", combine_blast, piden))
        #extract_top_n_hits(combine_blast + ".raw", top_num=max_hits, output=combine_blast, threads=threads, iden=iden)
    sh("cut -f1,2 {0} > {1}.raw".format(combine_blast, combine_ortho))
    uniq_ortho(combine_ortho + '.raw', combine_ortho)

    if runKs == 'T':
        if prefix2 == prefix1:
            combine_pep = prefix1 + '.pep'
            combine_cds = prefix1 + '.cds'
        else:
            combine_cds = combine_prefix + '.cds'
            combine_pep = combine_prefix + '.pep'
            sh("cat {0}.cds {1}.cds > {2}.cds".format(prefix1, prefix2, combine_prefix))
            sh("cat {0}.pep {1}.pep > {2}.pep".format(prefix1, prefix2, combine_prefix))
        kaks(combine_ortho, combine_cds, combine_pep, threads=threads, use_grid=use_grid, wait='F', model=model)
    return 0


def uniq_ortho(raw_ortho=None, uniq_ortho=None):
    """
    A simple script to find unique ortho regardless of ref and qry position.
    Two positional Args:
        raw_ortho: raw ortho file, read in
                   eg: At001g004t\tAt001g005t
                       At001g005t\tAt001g004t
        uniq_ortho: orthologs without redundant, write to
                   eg: At001g004t\tAt001g005t (last row in former eg is removed)
    Returns:
    """
    hit_table = {}
    with open(raw_ortho) as fh, open(uniq_ortho, 'w') as fh_out:
        for line in fh:
            mylist = line.rstrip().split()
            qry = mylist[0]
            ref = mylist[1]
            # Remove self ortholog like geneA-geneA, which makes no sense
            if qry == ref:
                continue
            tag1 = qry + "\t" + ref
            tag2 = ref + "\t" + qry
            if tag1 in hit_table or tag2 in hit_table:
                # Already recorded, continue
                continue
            else:
                hit_table[tag2] = 1
                hit_table[tag1] = 1
                # write tag1 is enough since tag1 and tag2 are reciprocal.
                fh_out.write(tag1 + "\n")


def rbhks(prefix1=None, prefix2=None, threads=4, runKs='T',
             use_grid='T', eval='1e-5', max_hits=10, model='YN', iden=0.2, output=''):
    r"""
    The Walker(2017) pipeline for calculating ortholog Ks
    Prerequisites: blast, KaKs_Calculator, ParaAT.pl
    ⭐️️The rbh_ks wrapper, execution eg:
        bsub python -m iga.evolution.ortho rbhks Cercis_chinensis Cercis_chinensis
    : param prefix1: like Cercis_chinensis. Require Cercis_chinensis.pep and Cercis_chinensis.gff3 exists
    : param prefix2: Can be the same with prefix1, which will perform intra comparison only.
    : param model: YN[default] or NG used for Ks estimation in ParaAT.pl (modified to accept model arg input)
    : param iden: min identity for blastp result (default 0.2 ), cannot be modified
    : param max_hits: 10, cannot be modified.
    :return:
    """
    combine_prefix = prefix1 + '.' + prefix2
    if output == '':
        workdir = 'work.{}'.format(combine_prefix)
    else:
        workdir = output
    #workdir_sh.format(workdir)
    mkdir(workdir)
    sh("cp {}.pep  {}/".format(prefix1, workdir))
    sh("cp {}.cds  {}/".format(prefix1, workdir))
    sh("cp {}.pep  {}/".format(prefix2, workdir))
    sh("cp {}.cds  {}/".format(prefix2, workdir))
    os.chdir(workdir)

    # input
    combine_blast = combine_prefix + '.blast'
    # output
    combine_ortho = combine_prefix + '.ortho'

    #Protein identity
    piden=20

    # prepare gff3 and blast
    if (not op.exists(combine_blast)):
        if not op.exists(combine_blast + '.raw'):
            blastp(prefix1 + ".pep", prefix2 + ".pep", threads=threads, output=combine_blast + ".raw",
                   use_grid='F', eval=eval)
        sh("awk '$3>={2}' {0} | onlyten.pl > {1}".format(combine_blast + ".raw", combine_blast, piden))
        #extract_top_n_hits(combine_blast + ".raw", top_num=max_hits, output=combine_blast, threads=threads, iden=iden)
    ("cut -f1,2 {0} > {1}.raw".format(combine_blast, combine_ortho))
    extract_reciprocal_best_hits(bln=combine_blast, output=combine_ortho + '.raw')
    uniq_ortho(combine_ortho + '.raw', combine_ortho)

    if runKs == 'T':
        if prefix2 == prefix1:
            combine_pep = prefix1 + '.pep'
            combine_cds = prefix1 + '.cds'
        else:
            combine_cds = combine_prefix + '.cds'
            combine_pep = combine_prefix + '.pep'
            sh("cat {0}.cds {1}.cds > {2}.cds".format(prefix1, prefix2, combine_prefix))
            sh("cat {0}.pep {1}.pep > {2}.pep".format(prefix1, prefix2, combine_prefix))
        kaks(combine_ortho, combine_cds, combine_pep, threads=threads, use_grid='F', wait='F', model=model)
    return 0


if __name__ == "__main__":
    emain()
