"""
code used in A188 project
"""
from collections import defaultdict

import pandas as pd
import numpy as np

from iga.apps.base import emain, logger, qsub, get_prefix, conda_act

# 0 ref fasta
# 1 qry fasta
nucmer_sh = r"""
# Whole genome alignment. Any other alignment can also be used.
WORKDIR={0}.{1}.nucmer
mkdir -p $WORKDIR
cd $WORKDIR
ln -s ../{0}
ln -s ../{1}
export PATH=/lustre/home/liuhui/bin/mummer4/bin:$PATH
# nucmer --maxmatch -c 100 -b 500 -l 50 {0} {1} 
nucmer --batch 1 -c 100 -b 500 -l 50 {0} {1}
# Remove small and lower quality alignments
delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta     
# Convert alignment information to a .TSV format as required by SyRI
show-coords -THrd out.filtered.delta > out.filtered.coords      
"""

syri_sh = r"""
SYRI=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/syri
PLOTSR=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/plotsr
python3 $SYRI -c out.filtered.coords -d out.filtered.delta -r {0} -q {1}
python3 $PLOTSR syri.out {0} {1} -H 8 -W 5
"""


def nucmer(ref=None, qry=None, threads=3):
    cmd = nucmer_sh.format(ref, qry)
    prefix = get_prefix(ref)
    prefix += get_prefix(qry)
    if len(ref.split('.')) > 2:
        chr_id = ref.split('.')[-1]
        prefix += chr_id
    qsub(cmd, cpus=threads, name='nucmer.'+prefix, sub=False)


def merge_nucmer_result():
    pass


def syri(ref=None, qry=None, threads=4):
    cmd = nucmer_sh.format(ref, qry) + '\nconda activate syri\n' + syri_sh.format(ref, qry)
    prefix = get_prefix(ref)
    prefix += get_prefix(qry)
    if len(ref.split('.')) > 2:
        chr_id = ref.split('.')[-1]
        prefix += chr_id
    qsub(cmd, cpus=threads, name='syri.'+prefix)


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

    def get_size(self):
        return int(self.end) - int(self.start) + 1


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

    def __init__(self, input_file='', type='bedpe'):
        # Nested storage system
        # level1: [dict] chromosome
        # level2: [list] Loci
        self.bedpe_db = defaultdict(list)
        self.type = type
        if input_file != '':
            self.read(input_file)

    def read(self, input_file):
        """
        Read a syri file and store
         A188_2  1307    3091    -       -       B73_2   12548   14331   SYNAL1  SYN1    SYNAL   -
        :return:
        """
        with open(input_file) as fh:
            for line in fh:
                mylist = line.rstrip().split()
                if self.type == 'syri':
                    (left_chr, left_start, left_end, undef, undef,
                     right_chr, right_start, right_end, name, name1, align_type, undef) = mylist
                    left_strand = '+'
                    right_strand = '+'
                elif self.type == 'bedpe':
                    (left_chr, left_start, left_end, right_chr, right_start, right_end, name,
                     undef, left_strand, right_strand) = mylist[0:10]
                this_lp = LociPE(left_chr, left_start, left_end, left_strand,
                                 right_chr, right_start, right_end, right_strand, name)
                self.bedpe_db[left_chr].append(this_lp)

    def stat(self):
        """
        stat bedpe size
        :return:
        """
        size_list_left = []
        size_list_right = []
        size_left_ins = []
        size_right_ins = []
        size_unknown_left = []
        size_unknown_right = []
        threshold = 50
        for chr_id in self.bedpe_db:
            chr_lp = self.bedpe_db[chr_id]
            for i, lp in enumerate(chr_lp):
                left_size = lp.left.get_size()
                right_size = lp.right.get_size()
                size_list_left.append(left_size)
                size_list_right.append(right_size)
                if left_size < threshold < right_size:
                    size_right_ins.append(right_size)
                elif left_size > threshold > right_size:
                    size_left_ins.append(left_size)
                elif left_size > threshold and right_size > threshold:
                    size_unknown_left.append(left_size)
                    size_unknown_right.append(right_size)

        pd.set_option('display.float_format', lambda x: '%.0f' % x)

        header = ['Left', 'Right', "Left Insertion", "Right Insertion", "Left Mosaic", "Right Mosaic"]
        tables = [size_list_left, size_list_right, size_left_ins, size_right_ins, size_unknown_left, size_unknown_right]

        # df = pd.DataFrame(np.array(tables), columns=header)
        # print(df.describe())
        for i, v in enumerate(header):
            print(header[i])
            if len(tables[i]) < 1:
                continue
            df = pd.DataFrame(tables[i])
            print(df.describe())

    def get_mosaic(self, outtable=''):
        """
        Get mosaic regions of this bedpe file
        :return:
        """
        # logger.debug("chrid {}".format(self.bedpe_db.keys()))
        complement_db = BedPE()
        for chr_id in self.bedpe_db:
            # logger.debug("chrid {}".format(chr_id))
            chr_lp = self.bedpe_db[chr_id]
            for i, lp in enumerate(chr_lp):
                if i > 0:
                    new_lp = LociPE(chr_lp[i - 1].left.chr, chr_lp[i - 1].left.end + 1, chr_lp[i].left.start - 1,
                                    chr_lp[i - 1].left.strand,
                                    chr_lp[i - 1].right.chr, chr_lp[i - 1].right.end + 1, chr_lp[i].right.start - 1,
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


def stat_bed(bedpe_file=None):
    """
    stat bedpe file
    :param bedpe_file:
    :return:
    """
    bedpe = BedPE(bedpe_file)
    bedpe.stat()


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
