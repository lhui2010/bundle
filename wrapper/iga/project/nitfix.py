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


def group2paralogs(orthogroup=None, max_group_size=10, start_col=3):
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


def group2orthologs(orthogroup=None, max_group_size=18, outdir='ortholog_split', start_col=3, threads=1):
    """
    12 min to finish
    Split groupt to parlogs
    Args:
        orthogroup: Orthogroups.tsv produced by OrthoFinder2.5.1
        threads: now deprecated
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
                     "Dipteryx_odorata_Pap",
                     "Pterodon_emarginatus_Pap",
                     "Zollernia_splendens_Pap",
                     "Dipteryx_alata",
                     "Eperua_falcata"]
    interest_list = ["Castanospermum_australe"]
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
    # logging.info(orthodb)
    # exit(1)
    # Tired of writing to multiple threads
    # with Pool(threads) as p:
    #     p.map(os.system, genblast_cmd_list)
    species_pairs = list(species_pairs)
    logging.info(species_pairs)
    # print(species_pairs)
    for spair in species_pairs:
        # print(spair)
        qry_dict = orthodb[spair[0]]
        ref_dict = orthodb[spair[1]]
        # logging.info(qry_dict)
        # logging.info(ref_dict)
        pair_name = "{}.{}".format(spair[0], spair[1])
        # if(len(qry_list) != len(ref_list)):
        #     logging.error('unequal length between {} and {}'.format(spair[0], spair[1]))

        iter_list = []
        for orthogroup in qry_dict:
            if orthogroup not in ref_dict:
                # logging.info("pass" + orthogroup)
                continue
            # https://stackoverflow.com/questions/12935194/permutations-between-two-lists-of-unequal-length
            # logging.info(qry_dict[orthogroup])
            # logging.info(ref_dict[orthogroup])
            else:
                iter_list.append([qry_dict[orthogroup], ref_dict[orthogroup]])
        # with Pool(threads) as pool:
        #     iter_result = pool.starmap(itertools.product, iter_list)
            iter_result = itertools.product(qry_dict[orthogroup], ref_dict[orthogroup])
            logging.info(list(iter_result))
        for iter_once_run in iter_result:
            # for ortho_pair in iter_result:
            for ortho_pair in iter_once_run:
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
    Returns:

    """
    if suffix1 == "":
        ortho1_list = ortho1.split('.')
        (ortho1_genera, ortho1_sp) = ortho1_list[0].split('_')
        suffix1 = ortho1_genera[:2].title() + ortho1_sp[:3]
    if suffix2 == "":
        ortho2_list = ortho2.split('.')
        (ortho2_genera, ortho2_sp) = ortho2_list[0].split('_')
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


emain()
