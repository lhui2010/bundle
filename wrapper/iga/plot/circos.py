"""
Circos plot wrappers
"""
import logging
import coloredlogs
from iga.apps.base import emain

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

def anchors_to_segdups(anchors=None, gene_bed=None):
    r"""
    Convert the anchors identified by MCScanX with gene location information (bed) to
    :param anchors:
    :param gene_bed:
    :return:
    """
    #Read bed in to dict

    #Read anchors and output anchors with bed information


def circos(fai=None, gene_gff='', gene_bed='', repeat_gff='', ortholog=''):


if __name__ == "__main__":
    emain()
