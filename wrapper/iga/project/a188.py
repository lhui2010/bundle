"""
code used in A188 project
"""
import copy
import re
from collections import defaultdict
from statistics import mean

import pandas as pd
from parse import parse

from iga.annotation.gff import Loci, Bed, GFF
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
        self.bedpe_db_right_chr = defaultdict(list)
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
                if self.type == 'lastz':
                    if len(self.bedpe_db[left_chr]) >= 2 and \
                            self.bedpe_db[left_chr][-1].left.start > this_lp.left.start:
                        break
                    elif len(self.bedpe_db[left_chr]) < 1 or \
                            self.bedpe_db[left_chr][-1].left.end < this_lp.left.start and \
                            self.bedpe_db[left_chr][-1].right.end < this_lp.right.end:
                        self.bedpe_db[left_chr].append(this_lp)
                else:
                    self.bedpe_db[left_chr].append(this_lp)
                    self.bedpe_db_right_chr[right_chr].append(this_lp)
                    # logging.debug(left_chr)
                    # logging.debug(right_chr)
                    # exit()

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
        count = 0
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
                    name = chr_lp[i - 1].right.name
                    if name == '.':
                        name = "SYN" + str(count)
                    new_lp = LociPE(chr_lp[i - 1].left.chr, left_start, left_end,
                                    chr_lp[i - 1].left.strand,
                                    chr_lp[i - 1].right.chr, right_start, right_end,
                                    chr_lp[i - 1].right.strand, "NOT" + name)
                    complement_db.bedpe_db[chr_id].append(new_lp)
                    count += 1
        complement_db.write_to_table(outtable)

    def write_to_table(self, table='', left_only=False, right_only=False):
        """
        write bedpe object into a table
        :param table: output table name
        :parameter left_only: print left bed only
        :parameter right_only: print right bed only
        :return:
        """
        result = ''
        for k in sorted(self.bedpe_db.keys()):
            # logger.debug(k)
            for i in self.bedpe_db[k]:
                if left_only:
                    result += i.left.get_line()
                elif right_only:
                    result += i.right.get_line()
                else:
                    result += i.get_line()
        if table != '':
            with open(table, 'w') as fh:
                fh.write(result)
        else:
            print(result, end='')

    def exists(self, bedpe_loci, wobble=100, type='ll', output='l'):
        """
        Test whether bedpe_loci also exists in in self.bedpe_db, with flexible boundary defined as wobble,
        and comparison type defined with type (left to left or right to right)
        :param bedpe_loci: the bedpe class object
        :param wobble: int, the boundaries can be flexible with plus or minus wobble
        :param type: the alignment type, could be ll, lr, rl, rr.
        :return:
        """
        abbrev = {"l": "left", "r": "right"}
        result = ''
        # bedpe_loci = LociPE()
        # bed_loci = Loci()
        if type[0] == 'l':
            # qry = bedpe_loci.abbrev[type[0]]
            qry = bedpe_loci.left
        else:
            qry = bedpe_loci.right

        if type[1] == 'l':
            pe_dict = self.bedpe_db
        else:
            pe_dict = self.bedpe_db_right_chr
        # logging.debug(qry.get_line())
        # logging.debug(qry.chr)
        # logging.debug(pe_dict.keys())
        # exit()
        for ref_pe in pe_dict[qry.chr]:
            # logging.debug(ref_pe.get_line())
            # exit()
            if type[1] == 'l':
                ref = ref_pe.left
            else:
                ref = ref_pe.right
            if abs(ref.start - qry.start) <= wobble and \
                    abs(ref.end - qry.end) <= wobble:
                if 'l' in output:
                    result += bedpe_loci.get_line().rstrip() + "\t"
                if 'r' in output:
                    result += ref_pe.get_line() + "\t"
                result = result.rstrip() + '\n'
        return result

    def hotspot(self, bedpe_loci, wobble=100, type='ll', output='l'):
        """
        Test whether bedpe_loci also exists in in self.bedpe_db, with flexible boundary defined as wobble,
        and comparison type defined with type (left to left or right to right)
        :param bedpe_loci: the bedpe class object
        :param wobble: int, the boundaries can be flexible with plus or minus wobble
        :param type: the alignment type, could be ll, lr, rl, rr.
        :return:
        """
        abbrev = {"l": "left", "r": "right"}
        start = ''
        end = ''
        # bedpe_loci = LociPE()
        # bed_loci = Loci()
        if type[0] == 'l':
            # qry = bedpe_loci.abbrev[type[0]]
            qry = bedpe_loci.left
        else:
            qry = bedpe_loci.right

        if type[1] == 'l':
            pe_dict = self.bedpe_db
        else:
            pe_dict = self.bedpe_db_right_chr
        # logging.debug(qry.get_line())
        # logging.debug(qry.chr)
        # logging.debug(pe_dict.keys())
        # exit()
        for ref_pe in pe_dict[qry.chr]:
            # logging.debug(ref_pe.get_line())
            # exit()
            if type[1] == 'l':
                ref = ref_pe.left
            else:
                ref = ref_pe.right
            if abs(ref.start - qry.start) <= wobble:
                # logging.debug('start found')
                if 'l' in output:
                    start += bedpe_loci.get_line().rstrip() + "\t"
                if 'r' in output:
                    start += ref_pe.get_line() + "\t"
                start = start.rstrip() + "\n"
                # logging.debug(start)
            if abs(ref.end - qry.end) <= wobble:
                # logging.debug('end found')
                if 'l' in output:
                    end += bedpe_loci.get_line().rstrip() + "\t"
                if 'r' in output:
                    end += ref_pe.get_line() + "\t"
                end = end.rstrip() + "\n"
                # logging.debug(end)
        # logging.debug([start, end])
        return [start, end]


def bedpe_stat(bedpe_file=None, short='F'):
    """
    stat bedpe file
    :param bedpe_file:
    :return:
    """
    bedpe = BedPE(bedpe_file)
    bedpe.stat(short)


def bed_stat(bed_file=None, short='F'):
    """
    stat bed file
    :param bed_file:
    :return:
    """
    bed = Bed(bed_file)
    bed.stat(short)


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
    # logging.debug('abc')
    if syriout == 'T':
        sh("grep SYNAL {0} > {0}.synal".format(synal_file))
        synal_file += '.synal'
    bedpe = BedPE(synal_file, type='syri')
    bedpe.get_mosaic()


def synal_to_paf(synal_file=None):
    """
    :param synal_file:
    :return:
    """
    # Description
    # 1 	string 	Query sequence name
    # 2 	int 	Query sequence length
    # 3 	int 	Query start (0-based; BED-like; closed)
    # 4 	int 	Query end (0-based; BED-like; open)
    # 5 	char 	Relative strand: "+" or "-"
    # 6 	string 	Target sequence name
    # 7 	int 	Target sequence length
    # 8 	int 	Target start on original strand (0-based)
    # 9 	int 	Target end on original strand (0-based)
    # 10 	int 	Number of residue matches
    # 11 	int 	Alignment block length
    # 12 	int 	Mapping quality (0-255; 255 for missing)
    bedpe = BedPE(synal_file, type='syri')
    for chr in bedpe.bedpe_db:
        for lp in bedpe.bedpe_db[chr]:
            paf_list = [lp.left.chr, 1000000, lp.left.start, lp.left.end, '+',
                        lp.right.chr, 1000000, lp.right.start, lp.right.end,
                        min(lp.left.end - lp.left.start, lp.right.end - lp.right.start),
                        max(lp.left.end - lp.left.start, lp.right.end - lp.right.start),
                        0]
            for i in range(0, len(paf_list)):
                paf_list[i] = str(paf_list[i])
            print("\t".join(paf_list))
    return 0


def split_paf(paf_file=None, bed_file=None, bin_size=1000000, offset='T'):
    """
    split minimap paf file into 1M-window seperated, for convenient plot with gggenome
    :param paf_file: arranged in A188 B73, B73 Mo17 order
    :param bed_file: the normal gene bed file. also transforming bed_file into split gff
    :return:
    """
    # Store split result in a list
    window_list = ['']
    # Store delimeter in a dict, like dict['a']=[1000000, 2000000]
    boundary_dict = {}
    bin_size = int(bin_size)
    longest_chr = 100000000000
    max_end = bin_size * int(longest_chr / bin_size)
    last_left_chr = ''
    last_right_chr = ''
    with open(paf_file) as fh:
        line = fh.readline()
        boundary_dict[line.split()[0]] = []
        boundary_dict[line.split()[5]] = []
        last_left_chr = line.split()[0]
        last_right_chr = line.split()[5]
        for i in range(bin_size, max_end, bin_size):
            boundary_dict[line.split()[0]].append(i)
    # Window number
    window_id = 0
    # Used to switch the two side
    side_list = ['left', 'right']
    known_side = 'left'
    unknown_side = 'right'
    last_unknown_end = 0
    with open(paf_file) as fh:
        for line in fh:
            this_line = defaultdict(dict)
            line_list = line.split()
            this_line['left']['chr'] = line_list[0]
            this_line['left']['start'] = line_list[2]
            this_line['left']['end'] = line_list[3]
            this_line['right']['chr'] = line_list[5]
            this_line['right']['start'] = line_list[7]
            this_line['right']['end'] = line_list[8]
            # (this_line['left']['chr'], undef1, this_line['left']['start'], this_line['left']['end'], undef2,
            #  this_line['right']['chr'], undef3, this_line['right']['start'], this_line['right']['end'], undef,
            #  undef, undef) = line.split()
            chr_id = this_line[known_side]['chr']
            # logging.debug(line)
            # logging.debug(chr_id)
            # logging.debug(this_line[known_side]['end'])
            # logging.debug(boundary_dict[chr_id][window_id])
            flag = False
            if last_left_chr == line_list[0] and last_right_chr == line_list[5]:
                flag = True
            if flag:
                # logging.debug(line)
                # logging.debug(chr_id)
                # logging.debug(window_id)
                # TODO: in current version start - offset can be negative
                if int(this_line[known_side]['end']) <= boundary_dict[chr_id][window_id]:
                    # window_list[window_id] += line
                    last_unknown_end = this_line[unknown_side]['end']
                else:
                    boundary_dict[this_line[unknown_side]['chr']].append(int(last_unknown_end))
                    window_id += 1
                    if len(window_list) <= window_id:
                        window_list.append('')
            else:
                # First judge whether this is a continuing block
                # Then judge which side is known side
                # But first finish the unfinished work
                boundary_dict[this_line[unknown_side]['chr']].append(int(last_unknown_end))
                if this_line[known_side]['chr'] not in boundary_dict:
                    # from A188 Mo17 to B73 Mo17
                    (known_side, unknown_side) = (unknown_side, known_side)
                if this_line[known_side]['chr'] not in boundary_dict:
                    logging.debug("Error: No known window offset in both columns")
                    exit(1)
                boundary_dict[this_line[unknown_side]['chr']] = []
                window_id = 0
                # window_list[window_id] += line
            if offset == 'T':
                if window_id > 0:
                    loff_set = boundary_dict[line_list[0]][window_id - 1]
                    roff_set = boundary_dict[line_list[5]][window_id - 1]
                else:
                    loff_set = 0
                    roff_set = 0
                line_list[3] = str(int(line_list[3]) - int(loff_set))
                line_list[2] = str(int(line_list[2]) - int(loff_set))
                line_list[8] = str(int(line_list[8]) - int(roff_set))
                line_list[7] = str(int(line_list[7]) - int(roff_set))
                line = "\t".join(line_list).rstrip() + "\n"
            window_list[window_id] += line
            last_left_chr = line_list[0]
            last_right_chr = line_list[5]
        # After loop end, appending last window offset
        boundary_dict[this_line[unknown_side]['chr']].append(int(last_unknown_end))

    for wd in range(0, len(window_list)):
        out_paf = '{}.{}.paf'.format(paf_file, wd)
        out_bed = '{}.{}.bed'.format(paf_file, wd)
        out_gff = '{}.{}.gff'.format(paf_file, wd)
        out_fake_fa = '{}.{}.fa'.format(paf_file, wd)
        out_fai = '{}.{}.fa.fai'.format(paf_file, wd)
        with open(out_paf, 'w') as fh:
            fh.write(window_list[wd])
        with open(out_bed, 'w') as fh:
            buffer = ''
            for k in boundary_dict:
                if wd == 0:
                    start = 0
                else:
                    start = boundary_dict[k][wd - 1]
                # logging.debug(k)
                # logging.debug(wd)
                # logging.debug(boundary_dict[k])
                try:
                    buffer += "{}\t{}\t{}\n".format(k, start, boundary_dict[k][wd])
                except IndexError:
                    logging.debug([k, start, wd])
            fh.write(buffer)
        sh('bedtools intersect -a {} -b {} -wb |cut -f4,5,6,7,8,9 |sort -k1,1 -k2,2n |uniq >{} '.format(
            out_bed, bed_file, out_bed + 'ist'))

        if offset == 'T':
            intersect_bed = Bed(out_bed + 'ist')
            with open(out_fai, 'w') as fh:
                for k in boundary_dict:
                    if wd == 0:
                        start = 0
                    else:
                        start = boundary_dict[k][wd - 1]
                    intersect_bed.change_offset(k, start)
                    try:
                        chr_size = boundary_dict[k][wd] - start
                        fh.write('{0}\t{1}\t{2}\t{3}\n'.format(k, chr_size, start, boundary_dict[k][wd]))
                    except IndexError:
                        logging.debug([k, start, wd])
                intersect_bed.write(out_bed + 'ist')
            sh('touch {}'.format(out_fake_fa))
        bed_to_gff(out_bed + 'ist', out_gff)
        # debug
        # break
    return 0


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


def split_by_tag(table=None, column=None):
    """
    split table by column
    :param table: 0-based
    :param colum:
    :return:
    """
    table_dict = defaultdict(str)
    column = int(column)
    with open(table) as fh:
        for line in fh:
            mylist = line.strip().split()
            table_dict[mylist[column]] += line
    for k in table_dict.keys():
        new_name = table + '.' + k
        with open(new_name, 'w') as fh:
            fh.write(table_dict[k])


def bedpe_intersect(bed1=None, bed2=None):
    """
    intersection of two bed pe file, like A188.B73.mosaic  and B73.Mo17.mosaic, default in left to right order
    :param bed1: A188.B73.mosaic
    :param bed2: B73.Mo17.mosaic
    :return:
    """
    bed1_buf = BedPE(bed1)
    bed2_buf = BedPE(bed2)

    bed1_out = bed1 + '.cut'
    bed2_out = bed2 + '.cut'

    bed1_buf.write_to_table(table=bed1_out, right_only=True)
    bed2_buf.write_to_table(table=bed2_out, left_only=True)

    sh("bedtools intersect -a {} -b {} > {}".format(bed1_out, bed2_out, bed1_out + bed2_out))

    intesect = Bed(bed1_out + bed2_out)

    intersect_sum = intesect.sum_size()
    print(intersect_sum)
    return intersect_sum


def bed_size(bed=None):
    """
    return bed size
    :param bed:
    :return:
    """
    bed = Bed(bed)
    print(bed.sum_size())
    return bed.sum_size()


def intersect_bedpe(bed1=None, bed2=None, type='ll', output='lr', wobble=100):
    """
    Get mosaic intersections
    :param bed1:
    :param bed2:
    :param type: could be ll lr rl rr. defining which side to be compared
    :param output: could be 'lr, l, r',defining pring both beds, left bed or right bed only
    :param wobble: the number of nucleotide that can be wobbled for the boundaries of mosaic region
    :return:
    """
    bed1_obj = BedPE(bed1)
    bed2_obj = BedPE(bed2)

    # self.bedpe_db[left_chr].append(this_lp)
    result = ''

    for chr in bed1_obj.bedpe_db:
        for bedpe_loci in bed1_obj.bedpe_db[chr]:
            search_result = bed2_obj.exists(bedpe_loci, wobble=int(wobble), type=type, output=output)
            if search_result != '':
                result += search_result
    print(result, end='')
    return result


def breakpoint_hotspot(bed1=None, bed2=None, type='ll', output='lr', wobble=100):
    """
    Get breakpoint hotspot between two bedPE files
    :param bed1:
    :param bed2:
    :param type: could be ll lr rl rr. defining which side to be compared
    :param output: could be 'lr, l, r',defining pring both beds, left bed or right bed only
    :param wobble: the number of nucleotide that can be wobbled for the boundaries of mosaic region
    :return:
    """
    bed1_obj = BedPE(bed1)
    bed2_obj = BedPE(bed2)

    # self.bedpe_db[left_chr].append(this_lp)
    start_all = ''
    end_all = ''

    for chr in bed1_obj.bedpe_db:
        for bedpe_loci in bed1_obj.bedpe_db[chr]:
            [start, end] = bed2_obj.hotspot(bedpe_loci, wobble=int(wobble), type=type, output=output)
            start_all += start
            end_all += end
    with open(bed1 + bed2 + '.start', 'w') as fh:
        fh.write(start_all)
    with open(bed1 + bed2 + '.end', 'w') as fh:
        fh.write(end_all)
    return 0


def bed_to_gff(bed=None, output=''):
    """
    Transform bed to gff format
    :param bed:
    :return:
    """
    bed_obj = Bed(bed)
    gff_obj = GFF()
    for k in bed_obj.bed_list:
        # logging.debug(k.name)
        if k.strand == '.':
            k.strand = '+'
        gff_obj.append(chr=k.chr, start=k.start, end=k.end, strand=k.strand, name=k.name)
    if output == '':
        gff_obj.print_out()
    else:
        with open(output, 'w') as fh:
            fh.write(gff_obj.to_str())


def annotate_block(anchor_file=None, qbed='', sbed=''):
    """
    alias: Ortho to bedpe
    Add ID, qstart, qend, sstart, send information for each syntenic block
    :param anchor_file: .anchor file generated by jcvi.compara.ortholog
    :param qbed: query bed file
    :param sbed: subject bed file
    :return:
    """
    qbed_obj = Bed(qbed)
    sbed_obj = Bed(sbed)
    # The ID for syntenic block
    syn_id = -1
    # The dictionary storing content for each syntenic block
    syn_content = defaultdict(str)
    syn_start = {}
    syn_end = {}
    with open(anchor_file) as fh:
        for line in fh:
            if line.startswith('#'):
                # starting of a syntenic block
                syn_id += 1
            else:
                (qgene, sgene, score) = line.rstrip().split()
                qloci = qbed_obj.select_name(qgene, format='loci')
                sloci = sbed_obj.select_name(sgene, format='loci')
                syn_content[syn_id] += qloci.get_line().strip() + "\t" + sloci.get_line()
                if syn_id not in syn_start:
                    syn_start[syn_id] = [qloci, sloci]
                if True:
                    syn_end[syn_id] = [qloci, sloci]
    for i in range(0, syn_id + 1):
        header = '#'
        qry_chr = syn_start[i][0].chr
        qry_start = syn_start[i][0].start
        qry_end = syn_end[i][0].end
        sub_chr = syn_start[i][1].chr
        sub_start = syn_start[i][1].start
        sub_end = syn_end[i][1].end
        strand = (sub_end - sub_start) * (qry_end - qry_start)
        if strand < 0:
            strand = '-'
        else:
            strand = '+'
        block_name = "Block{0}".format(i)
        tmp = "{}\t" * 8
        tmp = tmp.rstrip()
        header += tmp.format(qry_chr, qry_start, qry_end, sub_chr, sub_start, sub_end,
                             block_name, strand)
        content = syn_content[i]
        print(header + "\n" + content, end='')


def join_adjacent_bed(bed=None):
    """
    %prog join_adjacent_bed manual_reviewed_subgenome.bed > manual_reviewed_subgenome.joined.bed
    Input:
    1	6275775	6999144	ALN1
    1	7257371	15121456	ALN1
    Output:
    1 6275775 6275775 ALN1
    :param bed:
    :return:
    """
    bed_read = Bed(bed)
    loci = copy.copy(bed_read.bed_list[0])
    for i in range(1, len(bed_read.bed_list)):
        loci_i = bed_read.bed_list[i]
        if loci.chr == loci_i.chr and loci.name == loci_i.name:
            loci.start = min(loci.start, loci_i.start)
            loci.end = max(loci.end, loci_i.end)
        else:
            print(loci.get_line(), end='')
            loci = copy.copy(loci_i)
    print(loci.get_line(), end='')


def calc_percent_bed_intersect(bedwo=None):
    """
    %s A188_B73.total.SYN.mosaic.B73.bed.leaf_atac.wo > A188_B73.total.SYN.mosaic.B73.bed.leaf_atac.wo.perc
    Get relative position of intersection on Feature A
    Input:
        B73_1   2097669 2104245 NOTSYN26        B73_1   2104191 2105241 genic   54
        B73_1   2566205 2629140 NOTSYN34        B73_1   2628589 2628839 proximal        250
        B73_1   2636499 2746464 NOTSYN35        B73_1   2635661 2636711 genic   212
        B73_1   2892522 2894589 NOTSYN39        B73_1   2894533 2894758 proximal        56
    Output:
         B73_1   2097669 2104245 NOTSYN26        B73_1   2104191 2105241 genic   54      99.18%  100.00%
         B73_1   2566205 2629140 NOTSYN34        B73_1   2628589 2628839 proximal        250     99.12%  99.52%
         B73_1   2636499 2746464 NOTSYN35        B73_1   2635661 2636711 genic   212     0.00%   0.19%
         B73_1   2892522 2894589 NOTSYN39        B73_1   2894533 2894758 proximal        56      97.29%  100.00%
         B73_1   2973072 2974460 NOTSYN41        B73_1   2974176 2974451 genic   275     79.54%  99.35%
    :param bedwo:
    :return:
    """
    min_mosaic_size = 1000
    with open(bedwo) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            region_size = int(mylist[2]) - int(mylist[1])
            if region_size < min_mosaic_size:
                continue
            offset = int(mylist[1])
            rel_start = int(mylist[5]) - offset
            rel_end = int(mylist[6]) - offset
            start_percent = max(0, rel_start / region_size)
            end_percent = min(1, rel_end / region_size)
            print(line.rstrip() + "\t{:.2%}\t{:.2%}".format(start_percent, end_percent))


def percent_to_range(percent_file=None):
    r"""
    pass [a,b] to a a+1..b like format, for hist use
    Input:
        97.18% 100.00%
    Output:
        97
        98
        99
        100
    :param percent_file:
    :return:
    """
    import re
    with open(percent_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            for i in range(0, len(mylist)):
                mylist[i] = re.sub(r'\..*', '', mylist[i])
                mylist[i] = int(mylist[i])
            for a in range(mylist[0], mylist[1] + 1):
                print(a)


def breakpoint_screen(depth=None, highcutoff=100, lowcutoff=5):
    """
    %s A.depth.gz > A.breakpoint.txt
    Read depth file, and identify waterfall site that had a significant coverage drop
    :param highcutoff: larger than that is normal
    :param lowcutoff: lower than that is waterfall regions
    :param depth:
    :return:
    """
    import gzip
    higher_cutoff = highcutoff
    low_cutoff = lowcutoff
    with gzip.open(depth) as fh:
        prev_depth = 0
        prev_loci = 0
        prev_chr = ''
        buffer = None
        start_flag = False
        for line in fh:
            (chr_id, loci, depth) = line.decode().rstrip().split()
            if int(loci) % 10000000 == 0:
                logging.debug("Running on chr {} loci {}".format(chr_id, loci))
            depth = int(depth)
            if chr_id != prev_chr:
                if buffer is not None:
                    print(buffer.get_line(), end='')
            elif prev_depth > higher_cutoff and depth < low_cutoff:
                start_flag = 1
                buffer = Loci(chr_id, loci, loci, '.', '.', '.')
            elif start_flag and depth < low_cutoff:
                buffer.end = loci
            elif start_flag and depth > low_cutoff:
                print(buffer.get_line(), end='')
                buffer = None
                start_flag = False
            prev_depth = depth
            prev_loci = loci
            prev_chr = chr_id


def breakpoint_screen2(bam=None, add_name='F'):
    """
    %s test.bam > test.bam.breakpoint.txt
    add_name: [T|F]
    :return:
    """
    import pysam
    samfile = pysam.AlignmentFile(bam, "r")
    reads_all = samfile.fetch()
    buf = defaultdict(int)
    for read in reads_all:
        ###pysam's coordinate [0-based, 0-based), like bed, so have the following modifications
        # if type(read.reference_id) == int:
        #     read.reference_id += 1
        read.reference_start += 1
        #pysam 0.8.3, Future changed to latest version with reference_name
        try:
            reference_name = read.reference_name
        except AttributeError:
            reference_name = samfile.references[read.reference_id]
        ###
        if read.cigar[0][0] == 4 or read.cigar[0][0] == 5:
            print_buff = "{}\t{}\t{}".format(reference_name, read.reference_start, "Start")
            if add_name == 'T':
                print_buff += "\t{}".format(read.qname)
            print(print_buff)
            buf["{}\t{}\t{}".format(reference_name, read.reference_start, "Head")] += 1
        if read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5:
            print_buff = "{}\t{}\t{}".format(reference_name, read.reference_end, "End")
            if add_name == 'T':
                print_buff += "\t{}".format(read.qname)
            print(print_buff)
            buf["{}\t{}\t{}".format(reference_name, read.reference_end, "Tail")] += 1
    with open(bam + '.summary', 'w') as fh:
        for i in buf:
            fh.write("{}\t{}\n".format(i, buf[i]))


def add_depth_to_mosaic(mosaic_bedpe=None, bkptsum_l=None,
                        bkptsum_r=None, depth_l='', depth_r=''):
    """
    All is 1-based
    :param mosaic_bedpe:
    :param bkptsum:
    :param depth:
    :return:
    """
    mbe = BedPE(mosaic_bedpe)
    breakpoint_coverage_cutoff = 5
    bkptdb = defaultdict(int)
    bkptdb_right = defaultdict(int)
    with open(bkptsum_l) as fh:
        #1       37      Head    1
        for line in fh:
            mylist = line.rstrip().split()
            (chrid, loci, croptype, coverage) = mylist
            if int(coverage) < breakpoint_coverage_cutoff:
                continue
            bkptdb["_".join([chrid, loci])] += int(coverage)
    with open(bkptsum_r) as fh:
        #1       37      Head    1
        for line in fh:
            mylist = line.rstrip().split()
            try:
                (chrid, loci, croptype, coverage) = mylist
            except ValueError:
                logging.debug(line)
            if int(coverage) < breakpoint_coverage_cutoff:
                continue
            bkptdb_right["_".join([chrid, loci])] += int(coverage)
    breakpoint_flanking = 10
    for lchr in mbe.bedpe_db:
        logging.debug("left chromosome is {}".format(lchr))
        for lpe in mbe.bedpe_db[lchr]:
            lchr = re.sub(r'.*_', '', lpe.left.chr)
            rchr = re.sub(r'.*_', '', lpe.right.chr)
            # logging.debug('searching {}'.format(lpe.get_line()))
            etc = ''
            for left in (lpe.left.start, lpe.left.end):
                sum = 0
                for i in range(left - breakpoint_flanking, left + breakpoint_flanking):
                    try:
                        # logging.debug(bkptdb["_".join([lchr, i])])
                        # logging.debug(sum)
                        sum += bkptdb["_".join([lchr, str(i)])]
                    except KeyError:
                        pass
                etc += str(sum) + "\t"
            for right in (lpe.right.start, lpe.right.end):
                sum = 0
                for i in range(right - breakpoint_flanking, right + breakpoint_flanking):
                    try:
                        sum += bkptdb_right["_".join([rchr, str(i)])]
                    except KeyError:
                        pass
                etc += str(sum) + "\t"
            if etc != '':
                print(lpe.get_line().rstrip() + "\t" + etc.rstrip())

    # import gzip
    # with gzip.open(depth_l) as fh:
    #     for line in fh:
    #         (chr_id, loci, depth) = line.decode().rstrip().split()


if __name__ == "__main__":
    emain()
