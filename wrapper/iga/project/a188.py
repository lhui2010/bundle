"""
code used in A188 project
"""
from collections import defaultdict

from iga.apps.base import emain, logger


class Loci:
    """
    Loci object which could also be looked as bed object
    """

    def __init__(self, chr, start, end, name, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand


class LociPE:
    """
    Two loci objects
    """

    def __init__(self, left_chr, left_start, left_end, left_strand,
                 right_chr, right_start, right_end, right_strand, name):
        self.left = Loci(left_chr, left_start, left_end, name, left_strand)
        self.right = Loci(right_chr, right_start, right_end, name, right_strand)

    def get_line(self):
        """
        return string of this object
        :return:
        """
        result = "\t".join([self.left.chr, self.left.start, self.left.end, self.right.chr, self.right.start,
                            self.right.end, self.right.name, '.', self.left.strand, self.right.strand]) + "\n"
        return result


class BedPE:
    """
    parsing lastz or Syri result into BedPE object
    """

    def __init__(self, input_file='', type='syri'):
        # Nested storage system
        # level1: [dict] chromosome
        # level2: [list] Loci
        self.bedpe_db = defaultdict(list)
        if input_file != '':
            self.read(input_file, type)

    def read(self, input_file, type='syri'):
        """
        Read a syri file and store
         A188_2  1307    3091    -       -       B73_2   12548   14331   SYNAL1  SYN1    SYNAL   -
        :return:
        """
        with open(input_file) as fh:
            for line in fh:
                mylist = line.rstrip().split()
                if (type == 'syri'):
                    (left_chr, left_start, left_end, undef, undef,
                     right_chr, right_start, right_end, name, name1, align_type, undef) = mylist
                    left_strand = '+'
                    right_strand = '+'
                this_lp = LociPE(left_chr, left_start, left_end, left_strand,
                                 right_chr, right_start, right_end, right_strand, name)
                self.bedpe_db[left_chr].append(this_lp)

    def get_mosaic(self, outtable=''):
        """
        Get mosaic regions of this bedpe file
        :return:
        """
        logger.debug("chrid {}".format(self.bedpe_db.keys()))
        complement_db = BedPE()
        for chr_id in self.bedpe_db:
            logger.debug("chrid {}".format(chr_id))
            chr_lp = self.bedpe_db[chr_id]
            for i, lp in enumerate(chr_lp):
                if i > 0:
                    new_lp = LociPE(chr_lp[i - 1].left.chr, chr_lp[i - 1].left.start, chr_lp[i - 1].left.end,
                                    chr_lp[i - 1].left.strand,
                                    chr_lp[i - 1].right.chr, chr_lp[i - 1].right.start, chr_lp[i - 1].right.end,
                                    chr_lp[i - 1].right.strand, "NOT" + chr_lp[i - 1].right.name)
                    complement_db.bedpe_db[chr_id].append(new_lp)
        complement_db.write_to_table(outtable)

    def write_to_table(self, table=''):
        """
        write bedpe object into a table
        :param table:
        :return:
        """
        result = ''
        for k in sorted(self.bedpe_db.keys()):
            logger.debug(k)
            for i in self.bedpe_db[k]:
                result += i.get_line()
        if (table != ''):
            with open(table, 'w') as fh:
                fh.write(result)
        else:
            print(result, end='')


def synal_to_mosaic(synal_file=None):
    """
    %s syri.synal.txt > syri.unsynal.txt
    convert lastz result to bedpe like result by complementing syntenic regions
    Input eg:
    A188_2  1307    3091    -       -       B73_2   12548   14331   SYNAL1  SYN1    SYNAL   -
    A188_2  12474   16885   -       -       B73_2   31212   35640   SYNAL2  SYN2    SYNAL   -
    A188_2  16893   17371   -       -       B73_2   37259   37738   SYNAL3  SYN2    SYNAL   -
    A188_2  18572   20443   -       -       B73_2   39049   41277   SYNAL4  SYN2    SYNAL   -
    A188_2  25932   26692   -       -       B73_2   41807   42567   SYNAL5  SYN2    SYNAL   -
    Output eg:
    A188_2  3092    12473    -       -       B73_2   14332   31211   NONSYNAL1  SYN1    SYNAL   -
    A188_2  16886   16892   -       -       B73_2   35641   37258    NONSYNAL2  SYN2    SYNAL   -
    :param lastz:
    :return:
    """
    bedpe = BedPE(synal_file, type='syri')
    bedpe.get_mosaic()


if __name__ == "__main__":
    emain()
