"""
Circos plot wrappers
"""
import logging
import re

import coloredlogs

from iga.annotation.gff import BED
from iga.apps.base import emain, mkdir, Config, sh

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


def anchors_to_segdups(anchors=None, gene_bed=None):
    r"""
    Convert the anchors identified by MCScanX with gene location information (bed) to seg_dups used in circos
    :param anchors: anchor file built by jcvi, or the pair of orthologs described in the first two column
    :param gene_bed: gene bed files, could be generated by python -m iga.annotation.gff gff2bed gff [output: gff.bed]
    :return: the name of the segdups file, default is anchors + "seg_dups.txt"
    """
    #Read bed in to dict
    gene_bed_data = BED(gene_bed)
    #Read anchors and output anchors with bed information
    seg_dups_file = anchors + 'seg_dups.txt'
    seg_dups_buff = ''
    with open(anchors) as fh:
        for line in fh:
            if line.startswith('#'):
                seg_dups_buff += line.rstrip()
                continue
            mylist = line.split()
            (orthoA, orthoB) = mylist[:2]
            orthoA_bed = gene_bed_data.select_name(orthoA)
            orthoB_bed = gene_bed_data.select_name(orthoB)
            seg_line_list = [orthoA_bed.chr, str(orthoA_bed.start), str(orthoA_bed.end),
                orthoB_bed.chr, str(orthoB_bed.start), str(orthoB_bed.end)]
            seg_dups_buff += "  ".join(seg_line_list) + "\n"
    with open(seg_dups_file, 'w') as fh:
        fh.write(seg_dups_buff)
    return seg_dups_file


class Circos:
    """
    The circos class, support:
        Add plot structure like ggplot +
        Plot with available data
    Note:
        Because it's a little complicated, I used class instead of simple sh.format().
    """

    def __init__(self, karyotype_file):
        self.karyotype = ''
        self.link = ''
        self.ideogram = ''
        self.hist = []
        self.cytogenetic_bands = ''

        self.image_conf = Config('circos_image_conf')
        self.background_conf = Config('circos_background_conf')
        self.ideogram_conf = Config('circos_ideogram_conf')
        self.ticks_conf = Config('circos_ticks_conf')
        self.circos_conf = Config('circos')
        self.circos_conf.update('karyotype={}'.format(karyotype_file))

    def prepare_conf(self, hist_file=[], link_file=''):
        """
        prepare circos directory
        :return:
        """
        mkdir('data')
        mkdir('etc')
        self.image_conf.write_to_file('etc/image.conf')
        self.background_conf.write_to_file('etc/background.conf')
        self.ideogram_conf.write_to_file('ideogram.conf')
        self.ticks_conf.write_to_file('ticks.conf')
        self.circos_conf.write_to_file('circos.conf')

    def add_hist(self, hist_file='', plot_type=''):
        self.circos_conf.update('plots.plot')

    def plot(self):
        self.prepare_conf()
        sh('circos')


def fai_to_karyotype(fai=None):
    """
    fai to karyotype
    Input:
        chr01	14769264	2894042	70	71
    Output:
        chr - chr01 1 0 43445982 chr01
        chr - chr02 2 0 36083903 chr02
    :param fai:
    :return:
    """
    ky = ''
    ky_file = fai + ".karyotype"
    with open(fai) as fh:
        for line in fh:
            (chr_id, chr_size, undef, undef, undef) = line.split()
            try:
                chr_number = re.search(r'\d+', chr_id)[0]
                chr_number = re.sub(r'^0', '', chr_number)
            except TypeError:
                logger.error("Skipping non-chromosome contig {}".format(chr_id))
                continue
            ky += "chr - {0} {1} 0 {2} {0}\n".format(chr_id, chr_number, chr_size)
    with open(ky_file, 'w') as fh:
        fh.write(ky)
    return ky_file


def circos(fai=None, gene_gff='', gene_bed='', repeat_gff='', ortholog=''):
    karyotype_file = fai_to_karyotype(fai)
    circos_obj = Circos(karyotype_file)
    if ortholog != '':
        seg_dups_file = anchors_to_segdups(ortholog)
        #<<include links.conf>>
        #insert_link
    if gene_gff != '':
        gene_hist_file = gff_to_density(gene_gff)
    if repeat_gff != '':
        repeat_hist_file = gff_to_density(repeat_gff)
    circos_obj.plot()


if __name__ == "__main__":
    emain()
