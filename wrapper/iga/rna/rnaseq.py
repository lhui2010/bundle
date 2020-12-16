"""
rna-seq relevant utils
"""
from iga.apps.base import emain, get_prefix, logger
import pandas as pd
import numpy as np

from iga.plot.ggplot import pheatmap


def plot_exp_heatmap(table=None):
    """
    plot heatmap with merged table generated by merge_exp_table
    :param table:
    :return:
    """
    #First removing unrelated fields
    table2 = table + ".tmp"
    buff = ''
    with open(table) as fh:
        for line in fh:
            mylist = line.rstrip().split('\t')
            new_list = [mylist[0]] + mylist[6:]
            buff += "\t".join(new_list).rstrip() + "\n"
    with open(table2, 'w') as fh:
        fh.write(buff)
    pheatmap(table2, main="Expression Heatmap")
    return 0


def merge_exp_table(tables=None):
    """
    Merge expression table produced by stringtie into one with TPM
    :param tables:
    :return:
    """
    sample_list = []
    gene_dict = []
    z = pd.DataFrame()
    sub_list = ['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'TPM']
    for i, t in enumerate(tables):
        this_df = pd.read_table(t, sep='\t')
        logger.warning(this_df.columns)
        sub_df = this_df[list(sub_list)]
        t_prefix = get_prefix(t)
        new_sublist = sub_list.copy()
        new_sublist[new_sublist.index('TPM')] = t_prefix
        sub_df.columns = new_sublist
        if i == 0:
            z = sub_df
        else:
            z = z.merge(sub_df, left_on=sub_list[:6],
                    right_on=sub_list[:6], how='outer')
        #sample_list.append(get_prefix(t))
    z.to_csv('merged_expression.txt', sep="\t", index=False, na_rep='0.0')
    return z


if __name__ == "__main__":
    emain()
