"""
rna-seq relevant utils
"""
from iga.apps.base import emain, get_prefix
import pandas as pd
import numpy as np


def merge_expression_table(tables):
    sample_list = []
    gene_dict = []
    for t in tables:
        sample_list.append(get_prefix(t))
        pd.read_table(t)


if __name__ == "__main__":
    emain()
