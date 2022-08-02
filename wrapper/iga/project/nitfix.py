"""
scripts used in nitfix project
"""
import logging
import re
import coloredlogs
import ete3

from iga.annotation.gff import Bed
from iga.apps.base import emain, mkdir, sh, goto_workdir
import pandas as pd
import itertools
from collections import defaultdict
import os.path as op

from multiprocessing import Pool

import os
from ete3 import Tree

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

symbiosis_gene_fam_sh = r"""

RCG<-read.table("Orthogroups.GeneCount.symbiosis.tsv.count", header = T, row.names = 1, sep="\t")

library("pheatmap")
library("RColorBrewer")

# mycol =colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)

RCGm <- as.matrix(RCG)

RCG <- ifelse(RCGm>5, 5, RCGm)

mycol = brewer.pal(n=6, name="PuBu")

#
# col_size = length(colnames(RCG))
#
# RCG <- RCG[rowSums(RCG[,c(1:col_size)]) > 5]
#
# RCG <- log2(RCG/rowMeans(RCG))

prefix="Orthogroups.GeneCount.symbiosis.tsv.count"

pdf(paste(prefix, "pdf", sep='.'), w=10, h=30)

a = pheatmap(RCG, show_rownames=T, col = mycol,
       main = "Orthogroups.GeneCount.symbiosis.tsv.count", cluster_rows=T, cluster_cols=F,
       fontsize_col = 20, angle_col ="45",border_color = 'white')

"""


def unselect(ID=None, table=None):
    """
    like grep -vf, but much faster
    Args:
        ID:
        table:

    Returns:

    """
    id_dict = dict()
    with open(ID) as fh:
        for line in fh:
            id_dict[line.rstrip()] = 1
    with open(table) as fh:
        line = fh.readline()
        print(line.rstrip())
        for line in fh:
            if line.startswith('Sequence'):
                continue
            mylist = line.split()
            first_field = mylist[0].replace(">q_", "")
            (qry, ref) = first_field.split("_t_")
            if qry in id_dict or ref in id_dict:
                continue
            print(line.rstrip())


def format_ks_plot(ks_plot_output=None):
    """
    input: Umtiza_listeriana_Cae.Umtiza_listeriana_Cae.ortho.kaks.raw -> plot_ks_Umtiza_listeriana_Cae/paralogs-t300-mYN.kaks
    """
    fname_list = ks_plot_output.split('.')
    output_kaks = ks_plot_output.replace('.raw', '')
    output_ortho = ks_plot_output.replace('.kaks.raw', '')
    sp_name = fname_list[0]
    (genus, epithet) = sp_name.split('_')[:2]
    suffix = genus[:2] + epithet[:3]
    suffix = "_" + suffix.title()
    ortho_buff = ""
    ks_buff = ""
    with open(ks_plot_output) as fh:
        header = fh.readline()
        ks_buff += header
        for line in fh:
            if line.startswith('Sequence'):
                continue
            line = re.sub('>q_', '', line)
            line = re.sub('_t_', suffix + "-", line)
            line = re.sub('\t', suffix + "\t", line, 1)
            ks_buff += line
            # generate ortho file
            fields = line.split()[0].split('-')
            # logging.error(fields)
            # exit(1)
            ortho_buff += "{}\t{}\n".format(fields[0], fields[1])
    with open(output_kaks, 'w') as fh:
        fh.write(ks_buff)
    with open(output_ortho, 'w') as fh:
        fh.write(ortho_buff)


def group2paralogs(orthogroup=None, max_group_size=18, start_col=3):
    """
    12 min to finish
    Split groupt to parlogs
    Args:
        orthogroup: Orthogroups.tsv produced by OrthoFinder2.5.1
        start_col: 0-based indicating which column species started. 3 for N0.tsv and 1 for deprecated orthogroup.tsv
    Returns:
    """
    max_group_size = int(max_group_size)
    start_col = int(start_col)
    paralog_db = defaultdict(str)
    orthotable = pd.read_table(orthogroup, dtype=str)
    columns_len = len(orthotable.columns)
    for col in range(start_col, columns_len):
        species_name = orthotable.columns[col]
        this_species_groups = orthotable[species_name].to_list()
        for g in this_species_groups:
            try:
                ortho_gene_list = g.split(', ')
            except AttributeError:
                continue
            if (len(ortho_gene_list) <= 1 or len(ortho_gene_list) > max_group_size):
                continue
            else:
                ortho_pair_list = itertools.combinations(ortho_gene_list, 2)
                for ortho_pair in ortho_pair_list:
                    paralog_db[species_name] += ("\t".join(ortho_pair) + "\n")
    for species_name in paralog_db:
        with open("{}.{}.paralog".format(orthogroup, species_name), 'w') as fh:
            fh.write(paralog_db[species_name])


def group2orthologs(orthogroup=None, max_group_size=18, outdir='ortholog_split', start_col=3, threads=5, interest_sp=''):
    """
    12 min to finish
    Split groupt to parlogs
    Args:
        orthogroup: Orthogroups.tsv produced by OrthoFinder2.5.1
        threads: now OK, 3 or 5 is enough.
    Returns:
    """

    start_col = int(start_col)
    max_group_size = int(max_group_size)
    threads = int(threads)
    orthotable = pd.read_table(orthogroup, dtype=str, sep='\t')

    # Defing variables with prexisting values
    # Orthodb:
    #   species name(str)
    #       orthogroup(str)
    #           Ortholog of this group(list)
    orthodb = defaultdict(dict)

    columns_len = len(orthotable.columns)
    result_db = defaultdict(str)

    # Step1. GEt combinations only necessary
    # https://stackoverflow.com/questions/12935194/permutations-between-two-lists-of-unequal-length
    # cross comparison
    mkdir(outdir)
    species_pairs_raw = itertools.combinations(orthotable.columns[start_col:], 2)
    species_pairs = []
    interest_list = ["Andira_inermis_Pap",
                     "Dialium_schlechtneri_Dia",
                     "Goniorrhachis_marginata_Det",
                     "Umtiza_listeriana_Cae",
                     "Angylocalyx_braunii_Pap",
                     "Zollernia_splendens_Pap",
                     "Dipteryx_alata",
                     "Eperua_falcata"]
    if interest_sp != '':
        interest_list = [interest_sp]
    # interest_list2 = ["Castanospermum_australe"]
    for pair in species_pairs_raw:
        flag = 0
        for k in interest_list:
            if k in pair:
                flag = 1
        if (flag):
            species_pairs.append(pair)
    # species_pairs = itertools.product(interest_list2, interest_list)
    # logging.info(species_pairs)
    logging.info("Prepare of species pair complete, now writing to files...")

    for col in range(start_col, columns_len):
        species_name = orthotable.columns[col]
        this_species_groups = orthotable[species_name].to_list()
        for g in range(0, len(this_species_groups)):
            try:
                ortho_gene_list = this_species_groups[g].split(', ')
            except AttributeError:
                continue
            if (len(ortho_gene_list) <= 1 or len(ortho_gene_list) > max_group_size):
                continue
            else:
                orthodb[species_name][g] = ortho_gene_list
    species_pairs = list(species_pairs)
    logging.info(species_pairs)
    # print(species_pairs)
    for spair in species_pairs:
        # print(spair)
        qry_dict = orthodb[spair[0]]
        ref_dict = orthodb[spair[1]]
        pair_name = "{}.{}".format(spair[0], spair[1])
        # if(len(qry_list) != len(ref_list)):
        #     logging.error('unequal length between {} and {}'.format(spair[0], spair[1]))

        iter_list = []
        for orthogroup in qry_dict:
            if orthogroup not in ref_dict:
                # logging.info("pass" + orthogroup)
                continue
            # https://stackoverflow.com/questions/12935194/permutations-between-two-lists-of-unequal-length
            else:
                iter_list.append([qry_dict[orthogroup], ref_dict[orthogroup]])
        with Pool(threads) as pool:
            iter_result = pool.starmap(itertools.product, iter_list)
            # iter_result = itertools.product(qry_dict[orthogroup], ref_dict[orthogroup])
            # logging.info(list(iter_result))
        for iter_once_run in list(iter_result):
            for ortho_pair in iter_once_run:
            # for ortho_pair in iter_once_run:
            #     logging.info(ortho_pair)
                result_db[pair_name] += ("\t".join(ortho_pair) + "\n")
                # logging.info(pair_name)
        with open(op.join(outdir, pair_name + ".ortho"), 'w') as fh:
            fh.write(result_db[pair_name])
    # mkdir(outdir)
    # for pair_name in result_db:
    #     with open(op.join(outdir, pair_name + ".ortho"), 'w') as fh:
    #         fh.write(result_db[pair_name])


def remove_tandem_ortho(ortho=None, bed=None, gap=20):
    """
    Input:
        Ortho: AT01g01010\tAt02g02020
        Bed: bed file
    Output:
        Ortho with tandem duplicates removed (proximal genes with seperated by less than 20 genes)
        STDOUT!
    """
    bed_file = Bed(bed)
    # bed_file.add_rank()
    with open(ortho) as fh:
        for line in fh:
            (orthoA, orthoB) = line.rstrip().split()
            try:
                rankA = bed_file.select_name(orthoA).rank
                rankB = bed_file.select_name(orthoB).rank
            except KeyError:
                continue
            if (abs(rankA - rankB) > gap):
                print(line.rstrip())


# 0 orthologA-C
# 1 orthologB-C
# 2 Pepfile
# 3 threads
# 4 orthogroup type (mcscanx or orthofinder)
comWGD_tree_sh = """

touch {0} && rm {0}
touch {1} && rm {1}
touch {2} && rm {2}
ln -s ../{0}
ln -s ../{1}
ln -s ../{2}

BIN=select_ortholog_1_to_2.pl

if [ {4} == "orthofinder" ]
then
    BIN=select_ortholog_1_to_2.orthogroup.pl
fi

$BIN {0} > {0}.selection
$BIN {1} > {1}.selection

selectItem.pl -k {0}.selection {1}.selection  > {0}.{1}.group.txt

group2fasta.py {0}.{1}.group.txt {2} pep

touch tree
rm -rf tree

mkdir -p tree

cp */*pep tree/

pushd tree

touch run.sh
rm run.sh
for i in  *.pep
do 

    echo "mafft $i > $i.aln; iqtree2 -B 1000 -s $i.aln --prefix $i.aln -m LG+G" >> run.sh
done

cat run.sh | parallel -j {3} {{}}

popd

"""


def comWGD_tree(ortho1=None, ortho2=None, pep=None, suffix1='', suffix2='', threads=20, suffix_outgroup='Cechi',
                mode='mcscanx'):
    r"""
    Args:
        submit to execute

        ortho1: aa.cerchs.ortho
        ortho2: bb.cercis.ortho where cercis is outgroup
        pep:  aa.bb.pep  contain peptides of aa, bb, and cechi
        suffix1: Abrus_alba -> Abalb
        suffix2:
        suffix_outgroup: Cechi by default
        mode: [mcscanx|orthofinder], default mcscanx
    InputEg:
        ==> Castanospermum_australe.Cercis_chinensis.ortho <==
        CASAUS_g31283.t1_Caaus	CECHI00001252-t1_Cechi
        CASAUS_g31287.t2_Caaus	CECHI00001259-t1_Cechi
        CASAUS_g31303.t1_Caaus	CECHI00001266-t1_Cechi
        CASAUS_g31307.t1_Caaus	CECHI00001274-t1_Cechi
        CASAUS_g31317.t1_Caaus	CECHI00001279-t1_Cechi
        CASAUS_g31332.t2_Caaus	CECHI00001295-t1_Cechi
        CASAUS_g31339.t1_Caaus	CECHI00001302-t1_Cechi
        CASAUS_g31342.t1_Caaus	CECHI00001307-t1_Cechi
        CASAUS_g31347.t1_Caaus	CECHI00001309-t1_Cechi
        CASAUS_g31348.t1_Caaus	CECHI00001311-t1_Cechi

        ==> Phaseolus_lunatus.Cercis_chinensis.ortho <==
        Pl01G0000288100.1.v1_Phlun	CECHI00012229-t1_Cechi
        Pl01G0000288200.1.v1_Phlun	CECHI00012232-t1_Cechi
        Pl01G0000288500.1.v1_Phlun	CECHI00012236-t1_Cechi
        Pl01G0000288600.3.v1_Phlun	CECHI00012237-t3_Cechi
        Pl01G0000288700.1.v1_Phlun	CECHI00012239-t1_Cechi
        Pl01G0000288800.1.v1_Phlun	CECHI00012241-t1_Cechi
        Pl01G0000288900.1.v1_Phlun	CECHI00012242-t1_Cechi
        Pl01G0000289000.1.v1_Phlun	CECHI00012243-t1_Cechi
        Pl01G0000289200.1.v1_Phlun	CECHI00012244-t2_Cechi
    bash_eg:
        for QRY in Castanospermum_australe.ortho Cladrastis_platycarpa.ortho
            do
                for REF in *.ortho
                do
                    if [ $QRY == $REF ]
                    then
                        continue
                    fi
                    cat /ds3200_1/users_root/yitingshuang/lh/fasta0613/${QRY％.ortho}.pep \
                        /ds3200_1/users_root/yitingshuang/lh/fasta0613/${REF％.ortho}.pep \
                        /ds3200_1/users_root/yitingshuang/lh/fasta0613/Cercis_chinensis.pep  > ${QRY}.${REF}.pep

                    bsub -R "span[hosts=1]" -q Q104C512G_X4  -o ${QRY}.${REF}.log -e ${QRY}.${REF}.err -n 10
                    "python -m iga.project.nitfix comWGD_tree ${QRY} ${REF} ${QRY}.${REF}.pep --threads 10"
                 done
             done
    # collect results
        cp workdir_count_tree_*/*WGD results/
        for i in Nissolia_schottii Castanospermum_australe Cladrastis_platycarpa Pisum_sativum Styphnolobium_japonicum;
        do count_WGD_gene_tree.pl ${i}.*WGD > ${i}.count;
        done
    Returns:

    """
    ortho1_list = ortho1.split('.')
    (ortho1_genera, ortho1_sp) = ortho1_list[0].split('_')
    ortho2_list = ortho2.split('.')
    (ortho2_genera, ortho2_sp) = ortho2_list[0].split('_')
    if suffix1 == "":
        suffix1 = ortho1_genera[:2].title() + ortho1_sp[:3]
    if suffix2 == "":
        suffix2 = ortho2_genera[:2].title() + ortho2_sp[:3]
    goto_workdir('count_tree', ortho1_list[0] + '.' + ortho2_list[0])

    cmd = comWGD_tree_sh.format(ortho1, ortho2, pep, threads, mode)

    sh(cmd)
    output_files = yanrui_count_tree(tree_dir='tree',
                                     suffix_outgroup=suffix_outgroup,
                                     suffixA=suffix1,
                                     suffixB=suffix2)
    logging.info("Output files: ")
    logging.info(output_files)


def yanrui_count_tree(tree_dir=None, suffix_outgroup='', suffixA='', suffixB=''):
    """
    输入一个文件，三个参数
    Count tree topology (C是外群):
    (C,((A,B),(A,B)))的有116个，共有WGD
    (C,((A,A),(B,B)))的有755个，分别WGD
    输出三个文件：
    结果会包含不同类型的树结构。
    不符合上述两种噢的输出到other里面
    Args:
        tree_dir:
        suffix_outgroup:
        prefix:

    Returns:
    """

    nwk_file = [i for i in os.listdir(tree_dir) if i.endswith('treefile')]

    outgroup = suffix_outgroup
    spA = suffixA
    spB = suffixB

    com_print = ""
    ind_print = ""
    other_print = ""
    for i in nwk_file:
        # path='tree/' + i
        path = op.join(tree_dir, i)
        try:
            t = Tree(path)
        except ete3.parser.newick.NewickError:
            logging.error(path)
            exit(1)
        t_leave_name = [ii for ii in t.iter_leaf_names()]
        out_group_name = [ii for ii in t_leave_name if ii.endswith(outgroup)][0]  ####找到外类群######
        t.set_outgroup(out_group_name)  ####设置外类群，改变系统发育树的拓扑结构#####
        list_ind = []
        list_com = []
        for k in t.traverse():
            if k.is_leaf():
                continue
            else:
                child = k.get_children()[0].name + k.get_children()[1].name  ####获取节点名字#####
                if spA in child and spB in child:
                    # (C,((A,B),(A,B))) type is common
                    list_com.append(child)
                elif spA in k.get_children()[0].name and spA in k.get_children()[1].name:
                    # (C,((A,A),(B,B))) type is independent
                    list_ind.append(child)
                elif spB in k.get_children()[0].name and spB in k.get_children()[1].name:
                    list_ind.append(child)
        if len(list_ind) == 2:
            ind_print += t.write() + "\n"
        elif len(list_com) == 2:
            com_print += t.write() + "\n"
        else:
            other_print += t.write() + "\n"

    output_ind = "{}.{}.{}.indWGD".format(outgroup, spA, spB)
    output_com = "{}.{}.{}.comWGD".format(outgroup, spA, spB)
    output_other = "{}.{}.{}.otherWGD".format(outgroup, spA, spB)
    with open(output_ind, "w") as fh_ind, \
            open(output_com, "w") as fh_com, \
            open(output_other, "w") as fh_other:
        fh_ind.write(ind_print)
        fh_com.write(com_print)
        fh_other.write(other_print)
    return [output_ind, output_com, output_other]

# def automatic_get_root
# import ete3
# from ete3 import Tree
#
#
# def get_out_group(gene_tree=None,suffix=None):
#     with open(suffix) as fh:
#         suffix_list = fh.readlines()
#     tt = Tree(gene_tree)
#
#     root = None
# # Only one root is needed
#     del_list = []
#
#     for s in suffix_list:
#         for k in t.traverse():
#             if k.is_leaf():
#                 if("_".s in k.name):
#                     if root is None:
#                         root = k
#                         logging.info("Rooting on: ", k.name)
#                         continue
#                     else:
#                         del_list
#

# def plotTD(orthoFile=None, all_bed=None, total_pep=None, flank=10):
#     """
#     Args:
#         orthoFile:
#             MtgeneA   Castanospermum.geneA
#             MtgeneB   Castanospermum.geneB
#         all_bed:
#             cat Mt.bed Cast.bed > total.bed
#         total_pep:
#             cat Mt.pep Cast.pep > total.pep
#     Calculate:
#         1. fasttree of MtGeneA, MtGeneB and all homologs.
#         2. find flanking upstream and downstream of ten genes of Mt gene A
#         3. Use gggenes to plot location of MtGeneA and MtGeneB
#     Returns:
#     """
#     # Step1. create rename dict and reformat genes file. Used in rename nwk, beds.
#     xx.rename
#     # Step2. generate tree plot and rename with rename dict
#     iqtree
#     # Step3. get flanking genes with bed tools.
#     bedtools
#     # Step4. Reformat bed by replace .*: and .*- in column1.
#     perl
#     # Step5. Rename bed with rename dict.


#0 gene list
#1 gene alias
# output:
# {1}.tre
correct_gene_age_sh1 = r"""
set -euxo pipefail
get_homo_ortho2.sh {0} > homo_ortho.txt
get_syn_ortho.sh {0} > syn_ortho.txt
cat homo_ortho.txt syn_ortho.txt > raw_ortho.txt
tree.sh raw_ortho.txt
ln -s raw_ortho.txt.fa.aln.tre {1}.tre
echo "Manual rooting for {1}.tre" 
"""

# Progressive rooting

# 0 tree file; -q is set to 0.1 by default
# Output: {0}.shrink.tre {0}.shrink.txt
tree_shrink_sh = r"""
set +e
source activate treeshrink
run_treeshrink.py  -o . -O {0}.shrink -q 0.1 -t {0}
sed -i "s/\s\+/\n/" {0}.shrink.txt
echo "{0}.shrink.txt"
source deactivate
"""

# 0 : genes need to be removed;
# 1: gene alias
# 2: threads
correct_gene_age_sh2 = r"""
#. is directory; root is suffix for tree; 2 is branch length; 1 is minimum taxa (not used currently); output is out dir
# cut_long_internal_branches.py . root 2 1 output > long_branch.txt
unselectItem.pl {0} raw_ortho.txt  >clean_ortho.txt
# tree_iq2.sh clean_ortho.txt
tree_raxml.sh clean_ortho.txt {2}
ln -s clean_ortho.txt.fa.aln.treefile {1}.clean.tre
"""

# 0 rooted tree; MtCLE36.clean.tre
# output
# {0}.dlcdp.locus.recon
# {0}.dlcdp.locus.tree
dlcpar_sh = r"""
source activate dlcpar
# wait 300s then kill dlcpar as it is extreme slow on complex situations.
timeout 300 bash $BD/bash_template/dlcpar.sh {0}
echo "{0}.dlcdp.locus.tree"
source deactivate
"""

# 0: recon file. VsENBP1-like.rooted.tre.dlcdp.locus.recon
# 1: locus tree file. VsENBP1-like.rooted.tre.dlcdp.locus.tree
# output:
# MedicagoGene  DuplicationHistory  Orthologues
# MetruABC      N5:MetruABC,MetruDD|MetruEFT;N12MetruABC|MetruDD    AtABC,OrABC.
check_overlap_sh = r"""
python xx.py {0} {1}
"""


""" deprecated
tt = ete3.Tree("VsENBP1-like.rooted.tre.dlcdp.locus.tree", format=1)

for t in tt.traverse():
    if t.name == "n4":
        break
        

tt.get_common_ancestor("mRNA_MtrunA17Chr1g0147021_Metru", "mRNA_MtrunA17Chr7g0275771_M
etru")



list(filter(lambda x:'Metru' in x, t.children[0].get_leaf_names()))
['mRNA_MtrunA17Chr1g0147021_Metru']

list(filter(lambda x:'Metru' in x, t.children[1].get_leaf_names()))
['mRNA_MtrunA17Chr7g0275771_Metru']



activate dlcpar 
dlcpar dp -s ${sp_tree} -S ${sp_map} ${gene_tree} --output_format 3t

activate treeshrink
$run_treeshrink.py -t VsENBP1-like.all.tre -q 0.1


"""


def correct_gene_age(gene=None, threads=20):
    """
    submit to execute
    eg:
    correct_gene_age geneA,geneB,geneC
    Args:
        gene:

    Returns:

    """
    if type(gene) is list:
        gene = " ".join(gene)
    elif ',' in gene:
        gene = gene.replace(',', ' ')
    gene_alias = os.path.basename(os.getcwd())

    # tt = Tree(raw_tree)
    outgroup_list = ["Orsat",
"Aqcoe",
"Vivin",
"Artha,Potri,Avcar",
"Daglo,Bemas,Casat,Paand",
"Myrub,Cavim",
"Potat",
"Cechi,Bavar,Lyrho,Sigla",
"Duorc",
"Zeins",
"Mipud,Faalb,Sesep,Chpum"]

    # debug = 1

    prev_step =0
    # Step1 Build initial fasttree
    raw_tree = gene_alias + '.tre'
    if not os.path.exists(raw_tree) or prev_step:
        cmd1 = correct_gene_age_sh1.format(gene, gene_alias)
        sh(cmd1)
        prev_step = 1

    # Step2 root fasttree
    root_tree = raw_tree + ".root"
    if not os.path.exists(root_tree) or prev_step:
        _progressive_root_tree(raw_tree, outgroup_list)
        prev_step = 1
    # if not debug:
    #     _progressive_root_tree(raw_tree, outgroup_list)
    # output raw_tree + ".root

    # Step3 cut longbranch
    longbranch_ids = root_tree + '.shrink.txt'
    if not os.path.exists(longbranch_ids) or prev_step:
        cmd2 = tree_shrink_sh.format(root_tree)
        sh(cmd2)
        prev_step = 1

    # Step4: new tree with iqtree2
    clean_tree = gene_alias + ".clean.tre"
    if not os.path.exists(clean_tree) or prev_step:
        logging.info("Running step 4")
        cmd3 = correct_gene_age_sh2.format(longbranch_ids, gene_alias, threads)
        sh(cmd3)
        prev_step = 1
    else:
        logging.info("Passing step 4")

    # Step5: root iqtree2
    root_clean_tree = clean_tree + ".root"
    if not os.path.exists(root_clean_tree) or prev_step:
        logging.debug(clean_tree)
        _progressive_root_tree(clean_tree, outgroup_list)
        prev_step = 1

    # Step5
    dlcpar_input = root_clean_tree + '.filter'
    if not os.path.exists(dlcpar_input) or prev_step:
        _pre_dlcpar(root_clean_tree)
        prev_step = 1

    # Step6
    # 0 rooted tree; MtCLE36.clean.tre
    # output
    # {0}.dlcdp.locus.recon
    # {0}.dlcdp.locus.tree
    recon_file = dlcpar_input + ".dlcdp.locus.recon"
    locus_tree = dlcpar_input + ".dlcdp.locus.tree"
    if not os.path.exists(recon_file) or prev_step:
        cmd4 = dlcpar_sh.format(dlcpar_input)
        sh(cmd4)

    # Step 7
    if not os.path.exists(gene_alias + ".dups") or prev_step:
        result = _get_dups(recon_file, locus_tree)
        with open(gene_alias + '.dups', 'w') as fh:
            fh.write(result)
    #print(result, end="")


    # rooted_tree = input()
    # cmd2 = correct_gene_age_sh2
    # sh(cmd2)


def _progressive_root_tree(tree_fn, outgroup_list):
    """
    Args:
        tree: ETE3 tree object
        outgroup_list: ["Orsat",

    Returns:

    """
    tree = Tree(tree_fn)
    tree_labels = tree.get_leaf_names()
    outgroup_name = []
    for og in outgroup_list:
        og_list = og.split(',')
        # split og into suffix
        for split_og in og_list:
            split_og_list = [ii for ii in tree_labels if ii.endswith(split_og)]  ####找到外类群######
            outgroup_name += split_og_list
        if len(outgroup_name) > 0:
            break
    if len(outgroup_name) > 0:
        if len(outgroup_name) > 1:
            mrca_node = tree.get_common_ancestor(outgroup_name)
            mrca_outgroup_descends = list(filter(lambda x: x.is_leaf(),
                                             mrca_node.get_descendants()))
            logging.info(mrca_outgroup_descends)
            logging.info(outgroup_name)
            if len(mrca_outgroup_descends) == len(outgroup_name):
                # is monophyly
                tree.set_outgroup(mrca_node)
                tree.write(format=1, outfile=tree_fn+".root")
                return 0
        elif len(outgroup_name) == 1:
            tree.set_outgroup(outgroup_name[0])
            # TODO: change format to 0 to allow output of support value
            tree.write(format=1, outfile=tree_fn + ".root")
            return 0
        else:
            # is polyphyly
            logging.error("Polytomy found in outgroup for outgroup {} in {}".format(og,tree_fn))
            exit(1)
            return 1
    # No outgroup found. Return 1
    logging.error("No outgroup found for {}".format(tree_fn))
    exit(1)
    return 1


def _pre_dlcpar(tree_file):
    """
    Remove paralogous branches (Larger than 10) that underwent GDs but do not contain Metru nodes.
    Args:
        tree_file:

    Returns:tree_file + .filter

    """
    cutoff = 10
    tree = Tree(tree_file)
    children_list1 = tree.get_children()
    flag = ""
    for s in children_list1:
        if s.is_leaf():
            continue
        else:
            children_list2 = s.get_children()
            for st in children_list2:
                leaf_list = st.get_leaf_names()
                leaf_name_str = ",".join(leaf_list)
                if "_Metru" not in leaf_name_str and len(leaf_list) > cutoff:
                    logging.debug("Removing node", leaf_name_str)
                    flag = leaf_name_str
                    # tree.remove_child(st)
                    # st.detach()
                    break
    # tree.write(format=1, outfile=tree_file + ".filter")
    outfile = tree_file + ".filter"
    sh('pxrmt -t {} -n "{}" > {}'.format(tree_file, flag, outfile))
    return 0


def _get_dups(recon_file, locus_tree):
    """
    Args:
        recon_file:    # {0}.dlcdp.locus.recon
        locus_tree:    # {0}.dlcdp.locus.tree
    Returns:
    """
    # logging.info(recon_file)
    # logging.info(locus_tree)
    tree = Tree(locus_tree, format=1)
    # Duplication nodes
    result = ""
    with open(recon_file) as fh:
        for line in fh:
            (gene_tree_node, sp_tree_node, type) = line.rstrip().split()
            if type == "dup" and re.search(r"N\d+", sp_tree_node):
                # duplication node and happens on internal
                this_node = tree&gene_tree_node
                # this_node = tree.search_nodes(name=sp_tree_node)[0]
                children0_metru = list(filter(lambda x: 'Metru' in x, this_node.children[0].get_leaf_names()))
                children1_metru = list(filter(lambda x: 'Metru' in x, this_node.children[1].get_leaf_names()))
                # logging.debug(line)
                # logging.debug(children0_metru)
                # logging.debug(children1_metru)
                if len(children0_metru) > 0 and len(children1_metru) > 0:
                    treecp = tree.copy()
                    cut_node = treecp&gene_tree_node
                    cut_node.detach()
                    outgroup_genes = treecp.get_leaf_names()
                    result += "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        sp_tree_node,
                        ",".join(children0_metru),
                        ",".join(children1_metru),
                        ",".join(this_node.children[0].get_leaf_names()),
                        ",".join(this_node.children[1].get_leaf_names()),
                        ",".join(outgroup_genes)
                    )
    if result == "" and len(list(filter(lambda x: 'Metru' in x, tree.get_leaf_names()))) == 1:
        # No duplication at all. Report all genes and the outgroup
        if tree.children[0].name.startswith("n"):
            outgroup = tree.children[1].name
        else:
            outgroup = tree.children[0].name
        result += "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        outgroup,
                        "-",
                        "-",
                        "-",
                        "-",
                        ",".join(tree.get_leaf_names()))

    return result


def get_dups(recon_file=None, locus_tree=None):
    """
    Input: dlcpar output
    Args:
        recon_file:    # {0}.dlcdp.locus.recon
        locus_tree:    # {0}.dlcdp.locus.tree
    Returns: print duplication result in STDOUT
    """
    result = _get_dups(recon_file, locus_tree)
    print(result)


if __name__ == "__main__":
    emain()
