"""
scripts used in nitfix project
"""
import logging
import re
import coloredlogs

from iga.annotation.gff import Bed
from iga.apps.base import emain, mkdir
import pandas as pd
import itertools
from collections import defaultdict
import os.path as op

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
            #generate ortho file
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
            if (len(ortho_gene_list) <=1 or len(ortho_gene_list) > max_group_size):
                continue
            else:
                ortho_pair_list = itertools.combinations(ortho_gene_list, 2)
                for ortho_pair in ortho_pair_list:
                    paralog_db[species_name] += ("\t".join(ortho_pair)+"\n")
    for species_name in paralog_db:
        with open("{}.{}.paralog".format(orthogroup, species_name), 'w') as fh:
            fh.write(paralog_db[species_name])


def group2orthologs(orthogroup=None, max_group_size=18, outdir='ortholog_split', start_col=3):
    """
    12 min to finish
    Split groupt to parlogs
    Args:
        orthogroup: Orthogroups.tsv produced by OrthoFinder2.5.1
    Returns:
    """
    start_col = int(start_col)
    max_group_size = int(max_group_size)
    ortho_db = defaultdict(str)
    orthotable = pd.read_table(orthogroup, dtype=str, sep='\t')
    columns_len = len(orthotable.columns)
    result_db = defaultdict(str)
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
                ortho_db[species_name] = ortho_gene_list
    #https://stackoverflow.com/questions/12935194/permutations-between-two-lists-of-unequal-length
    #cross comparison
    mkdir(outdir)
    species_pairs_raw = itertools.combinations(orthotable.columns[start_col:], 2)
    species_pairs = []
    interest_list = ["Andira_inermis_Pap", "Dialium_schlechtneri_Dia", "Goniorrhachis_marginata_Det", "Umtiza_listeriana_Cae", "Angylocalyx_braunii_Pap", "Dipteryx_odorata_Pap", "Pterodon_emarginatus_Pap", "Zollernia_splendens_Pap", "Dipteryx_alata", "Eperua_falcata"]

    for pair in species_pairs_raw:
        flag = 0
        for k in interest_list:
            if k in pair:
                flag = 1
        if(flag):
            species_pairs.append(pair)
    logging.info(species_pairs)
    logging.info("Prepare of species pair complete, now writing to files...")
    for spair in species_pairs:
        qry_list = orthotable[spair[0]]
        ref_list = orthotable[spair[1]]
        pair_name = "{}.{}".format(spair[0], spair[1])
        if(len(qry_list) != len(ref_list)):
            logging.error('unequal length between {} and {}'.format(spair[0], spair[1]))
        for row_id in range(0, len(qry_list)):
            #https://stackoverflow.com/questions/12935194/permutations-between-two-lists-of-unequal-length
            logging.info(orthotable[spair[0]][row_id])
            logging.info(orthotable[spair[1]][row_id])
            iter_result = itertools.product(orthotable[spair[0]][row_id], orthotable[spair[1]][row_id])
            for ortho_pair in iter_result:
                result_db[pair_name] += ("\t".join(ortho_pair)+"\n")
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
    #bed_file.add_rank()
    with open(ortho) as fh:
        for line in fh:
            (orthoA, orthoB) = line.rstrip().split()
            try:
                rankA = bed_file.select_name(orthoA).rank
                rankB = bed_file.select_name(orthoB).rank
            except KeyError:
                continue
            if(abs(rankA - rankB) > gap):
                print(line.rstrip())


emain()