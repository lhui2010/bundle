"""
scripts used in nitfix project
"""
import logging
import re
import coloredlogs
from iga.apps.base import emain
import pandas as pd
import itertools
from collections import defaultdict

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


def group2paralogs(orthogroup=None, max_group_size=10):
    """
    Split groupt to parlogs
    Args:
        orthogroup:
    Returns:
    """
    paralog_db = defaultdict(str)
    orthotable = pd.read_table(orthogroup)
    columns_len = len(orthotable.columns)
    for col in range(1, columns_len - 1):
        species_name = orthotable.columns[col]
        this_species_groups = orthotable[species_name].to_list()
        for g in this_species_groups:
            ortho_gene_list = g.split(', ')
            if (len(ortho_gene_list) <=1 or len(ortho_gene_list) > max_group_size):
                continue
            else:
                ortho_pair_list = itertools.combinations(ortho_gene_list, 2)
                for ortho_pair in ortho_pair_list:
                    paralog_db[species_name] += ("\t".join(ortho_pair)+"\n")
    for species_name in paralog_db:
        with open("{}.{}.paralog".format(orthogroup, species_name), 'w') as fh:
            fh.write(paralog_db[species_name])


emain()
