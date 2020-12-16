"""
rna-seq relevant utils
"""
from iga.apps.base import emain, get_prefix
import pandas as pd
import numpy as np


def merge_exp_table(tables=None):
    """
    Merge expression table produced by stringtie into one with TPM
    :param tables:
    :return:
    """
    sample_list = []
    gene_dict = []
    z = pd.DataFrame()
    sub_list = ('Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'TPM')
    for i, t in enumerate(tables):
        this_df = pd.read_table(t, sep='\t')
        this_df
        sub_df = this_df[sub_list]
        t_prefix = get_prefix(t)
        new_sublist = list(sub_list)
        new_sublist[new_sublist.index('TPM')] = t_prefix
        sub_df.columns = new_sublist
        if i == 0:
            z = sub_df
        else:
            z.merge(sub_df, left_on=['Gene', 'mRNA', 'chr', 'strand', 'start', 'end'],
                    right_on=['Gene', 'mRNA', 'chr', 'strand', 'start', 'end'], how='outer')
        #sample_list.append(get_prefix(t))
    z.to_csv('merged_expression.txt', sep="\t")


if __name__ == "__main__":
    emain()
