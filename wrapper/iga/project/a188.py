"""
code used in A188 project
"""
from collections import defaultdict
from statistics import mean

import pandas as pd
from parse import parse

from iga.annotation.gff import Loci
from iga.apps.base import emain, qsub, get_prefix, sh

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# # 0 ref fasta
# # 1 qry fasta
# nucmer_sh = r"""
# # Whole genome alignment. Any other alignment can also be used.
# WORKDIR={0}.{1}.nucmer
# mkdir -p $WORKDIR
# cd $WORKDIR
# ln -s ../{0}
# ln -s ../{1}
# export PATH=/lustre/home/liuhui/bin/mummer4/bin:$PATH
# # nucmer --maxmatch -c 100 -b 500 -l 50 {0} {1}
# nucmer --batch 1 -c 100 -b 500 -l 50 {0} {1}
# # Remove small and lower quality alignments
# delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta
# # Convert alignment information to a .TSV format as required by SyRI
# show-coords -THrd out.filtered.delta > out.filtered.coords
# """


# syri_sh = r"""
# SYRI=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/syri
# PLOTSR=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/plotsr
# python3 $SYRI -c out.filtered.coords -d out.filtered.delta -r {0} -q {1}
# python3 $PLOTSR syri.out {0} {1} -H 8 -W 5
# """


# def nucmer(ref=None, qry=None, threads=3):
#     cmd = nucmer_sh.format(ref, qry)
#     prefix = get_prefix(ref)
#     prefix += get_prefix(qry)
#     if len(ref.split('.')) > 2:
#         chr_id = ref.split('.')[-1]
#         prefix += chr_id
#     qsub(cmd, cpus=threads, name='nucmer.'+prefix, sub=False)


# def merge_nucmer_result():
#     pass


# def syri(ref=None, qry=None, threads=4):
#     cmd = nucmer_sh.format(ref, qry) + '\nconda activate syri\n' + syri_sh.format(ref, qry)
#     prefix = get_prefix(ref)
#     prefix += get_prefix(qry)
#     if len(ref.split('.')) > 2:
#         chr_id = ref.split('.')[-1]
#         prefix += chr_id
#     qsub(cmd, cpus=threads, name='syri.'+prefix)


syri_sh = r"""
REF={0}
QRY={1}
WORKDIR={0}.{1}.syri
mkdir -p $WORKDIR
cd $WORKDIR
ln -s ../{0}
ln -s ../{1}
minimap2 -ax asm5 --eqx $REF $QRY | samtools view -bS  > $REF.$QRY.bam
SYRI=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/syri
PLOTSR=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/plotsr

python3 $SYRI -c $REF.$QRY.bam -r $REF -q $QRY -k -F B --lf $REF.$QRY.log --prefix $REF.$QRY.
# python3 $SYRI -c out.filtered.coords -d out.filtered.delta -r {0} -q {1}

for i in .snps.txt .notAligned.txt .ctxOut.txt .sv.txt .synOut.txt .dupOut.txt .invDupOut.txt .invTLOut.txt .TLOut.txt .invOut.txt
do 
    rm $REF.$QRY$i
done
python3 $PLOTSR $REF.$QRY.syri.out {0} {1} -H 8 -W 5
"""


def syri(ref=None, qry=None, threads=6):
    cmd = 'conda activate syri' + syri_sh.format(ref, qry)
    prefix = get_prefix(ref)
    prefix += get_prefix(qry)
    if len(ref.split('.')) > 2:
        chr_id = ref.split('.')[-1]
        prefix += chr_id
    qsub(cmd, cpus=threads, name='syri.' + prefix)


class LociPE:
    """
    Two loci objects
    """

    def __init__(self, left_chr, left_start, left_end, left_strand,
                 right_chr, right_start, right_end, right_strand, name):
        self.left = Loci(left_chr, left_start, left_end, name, '.', left_strand)
        self.right = Loci(right_chr, right_start, right_end, name, '.', right_strand)

    def get_line(self):
        """
        return string of this object
        :return:
        """
        result = "\t".join([self.left.chr, str(self.left.start), str(self.left.end),
                            self.right.chr, str(self.right.start), str(self.right.end),
                            self.right.name, '.', self.left.strand, self.right.strand]) + "\n"
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
                if line.startswith('#'):
                    continue
                mylist = line.rstrip().split()
                if self.type == 'syri':
                    (left_chr, left_start, left_end, undef, undef,
                     right_chr, right_start, right_end, name, name1, align_type, undef) = mylist
                    left_strand = '+'
                    right_strand = '+'
                    # Default input for syri is 1-based for both start and end, so have a change to bedpe type
                    right_start = int(right_start) - 1
                    left_start = int(left_start) - 1
                elif self.type == 'bedpe':
                    # [0-based start, 0-based end)
                    (left_chr, left_start, left_end, right_chr, right_start, right_end, name,
                     undef, left_strand, right_strand) = mylist[0:10]
                elif self.type == 'lastz':
                    # lastz is also compatible with bed format: [0-based start, 0-based end)
                    # name1  zstart1 end1 name2   strand2 zstart2+  end2+ identity idPct coverage covPct  cigarx-
                    (left_chr, left_start, left_end, right_chr, right_strand, right_start, right_end,
                     identity, idcPct, coverage, covPct, cigarx) = mylist
                    left_strand = '.'
                    name = '.'

                this_lp = LociPE(left_chr, left_start, left_end, left_strand,
                                 right_chr, right_start, right_end, right_strand, name)
                self.bedpe_db[left_chr].append(this_lp)

    def stat(self, short):
        """
        stat bedpe size
        :param short:, whether to use short format, T is use
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
                left_size = lp.left.get_size(bed_format=True)
                right_size = lp.right.get_size(bed_format=True)
                if left_size == 1:
                    # Skip loci with size of 1
                    left_size = 0
                if right_size == 1:
                    right_size = 0
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

        header = ['Left', 'Right', "LeftIns", "RightIns", "Left Mosaic", "Right Mosaic"]
        tables = [size_list_left, size_list_right, size_left_ins, size_right_ins, size_unknown_left, size_unknown_right]

        # df = pd.DataFrame(np.array(tables), columns=header)
        # print(df.describe())
        if short == 'T':
            print("{}\t{}".format(sum(tables[0]), sum(tables[1])))
        else:
            print("\t", end='')
            for h in header:
                print("{:<15}".format(h), end='')
            print('')

            print("max", end="\t")
            for i in range(0, len(header)):
                print("{:<15}".format(max(tables[i])), end='\t')
            print('')

            print("mean", end="\t")
            for i in range(0, len(header)):
                print("{:<15}".format(int(mean(tables[i]))), end='\t')
            print('')

            print("sum", end="\t")
            for i in range(0, len(header)):
                print("{:<15}".format(sum(tables[i])), end='\t')
            print('')
        return 0

        # for i, v in enumerate(header):
        #     print(header[i])
        #     if len(tables[i]) < 1:
        #         continue
        #     df = pd.DataFrame(tables[i])
        #     print(df.describe())

    def get_mosaic(self, outtable=''):
        """
        Get mosaic regions of this bedpe file
        CHANGE0128: Use increment alignment only
        :return:
        """
        # logger.debug("chrid {}".format(self.bedpe_db.keys()))
        complement_db = BedPE()
        for chr_id in self.bedpe_db:
            # logger.debug("chrid {}".format(chr_id))
            chr_lp = self.bedpe_db[chr_id]
            for i, lp in enumerate(chr_lp):
                if i > 0:
                    left_start = chr_lp[i - 1].left.end + 1
                    left_end = chr_lp[i].left.start - 1
                    right_start = chr_lp[i - 1].right.end + 1
                    right_end = chr_lp[i].right.start - 1
                    if left_end <= left_start:
                        left_end = left_start
                    if right_end <= right_start:
                        right_end = right_start
                    # From now, 1-based is transformed into 0-based start and 1-based end.
                    new_lp = LociPE(chr_lp[i - 1].left.chr, left_start, left_end,
                                    chr_lp[i - 1].left.strand,
                                    chr_lp[i - 1].right.chr, right_start, right_end,
                                    chr_lp[i - 1].right.strand, "NOT" + chr_lp[i - 1].right.name)
                    # TODO just a tempory fix for nonincrement alignment
                    if len(complement_db.bedpe_db[chr_id]) > 2 and \
                            (complement_db.bedpe_db[chr_id][-1].left.start > new_lp.left.start or
                             complement_db.bedpe_db[chr_id][-1].right.start > new_lp.right.start) and \
                            (complement_db.bedpe_db[chr_id][-2].left.start < new_lp.left.start and
                             complement_db.bedpe_db[chr_id][-2].right.start < new_lp.right.start):
                        complement_db.bedpe_db.pop()
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
        if table != '':
            with open(table, 'w') as fh:
                fh.write(result)
        else:
            print(result, end='')


def stat_bed(bedpe_file=None, short='F'):
    """
    stat bedpe file
    :param bedpe_file:
    :return:
    """
    bedpe = BedPE(bedpe_file)
    bedpe.stat(short)


def synal_to_mosaic(synal_file=None, syriout='F'):
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
    :param synal_file: Input alignment file
    :param syriout: [F/T] whether this is syri.out file
    :return:
    """
    if syriout == 'T':
        sh("grep SYNAL {0} > {0}.synal".format(synal_file))
        synal_file += '.synal'
    bedpe = BedPE(synal_file, type='syri')
    bedpe.get_mosaic()


def lastz_to_mosaic(lastz_file=None):
    """
    %s lastz.txt(sort -k1,1) > lastz.mosaic.txt
    :param lastz_file: Input alignment file
    :return: STDOUT
    """
    bedpe = BedPE(lastz_file, type='lastz')
    bedpe.get_mosaic()


def mosaic_ratio(fai=None, stat=None):
    """
    Calculate mosaic ratio based on mosaic region size and chromosome size
    :param fai: chromosome size
    :param stat: mosaic ratio (format: left_mosaic_size\tright_mosaic_size)
    :return:
    """
    chr_size = {}
    with open(fai) as fh:
        for line in fh:
            mylist = line.split()
            chr_size[mylist[0]] = int(mylist[1])
    with open(stat) as fh:
        line = fh.readline()
        (left_mosaic_size, right_mosaic_size) = line.strip().split()
    # PH207.genome.chr10.W22.genome.chr10.syri.out.mosaic.stat
    mylist = parse("{}.genome.chr{}.{}.genome.{}.syri.out.mosaic.stat", stat)
    tag = "{}-{}-{}".format(mylist[0], mylist[2], mylist[1])
    left_chr = "{}_{}".format(mylist[0], mylist[1])
    right_chr = "{}_{}".format(mylist[2], mylist[1])
    ratio = (int(left_mosaic_size) + int(right_mosaic_size)) / (chr_size[left_chr] + chr_size[right_chr])
    print("{}\t{}\t{}\t{}\t{}\t{:.1%}".format(tag, left_mosaic_size, right_mosaic_size,
                                              chr_size[left_chr], chr_size[right_chr], ratio))


def chromosome_level_ratio(stat=None):
    """
    Run this function after
    `$for i in *stat; do python -m iga.project.a188 mosaic_ratio total.fai $i >> mosaic_ratio_stat.txt; done`
    Input eg:
        A188-B73-10     66760136        70522917        146705344       150982314       46.1%
        A188-Mo17-10    70611181        72223243        146705344       149041351       48.3%
        A188-PH207-10   76993877        74673644        146705344       145727586       51.9%
        A188-SK-10      73901986        75129013        146705344       148363080       50.5%
    :param stat: mosaic_ratio_stat.txt
    :return:
    """
    cmp_dict1 = defaultdict(int)
    cmp_dict2 = defaultdict(int)
    cmp_dict3 = defaultdict(int)
    cmp_dict4 = defaultdict(int)
    with open(stat) as fh:
        for line in fh:
            (tag, l_mosaic_size, r_mosaic_size, l_chr_size, r_chr_size, ratio) = line.split()
            (l_id, r_id, chr_id) = tag.split('-')
            chr_tag = l_id + '-' + r_id
            cmp_dict1[chr_tag] += int(l_mosaic_size)
            cmp_dict2[chr_tag] += int(r_mosaic_size)
            cmp_dict3[chr_tag] += int(l_chr_size)
            cmp_dict4[chr_tag] += int(r_chr_size)
    for k in cmp_dict1:
        ratio = (cmp_dict1[k] + cmp_dict2[k]) / (cmp_dict3[k] + cmp_dict4[k])
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2%}".format(k, cmp_dict1[k], cmp_dict2[k],
                                                        cmp_dict3[k], cmp_dict4[k], ratio))


if __name__ == "__main__":
    emain()
