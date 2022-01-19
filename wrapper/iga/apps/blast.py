"""
blast related wrappers
"""
from collections import defaultdict

from iga.apps.base import emain, bsub, waitjob, sh

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# 0 ref
# 1 qry
# 2 output
# 3 eval
# 4 threads
# 5 out
blastp_sh = r"""
makeblastdb -in {0} -dbtype prot
blastp  -num_threads {4} -db {0} -query {1}  -evalue {3} -outfmt 6 -out {5}
"""


def blastp(ref=None, qry=None, threads=5, eval=1e-5, use_grid='T', output='-'):
    r"""
    blastp wrapper, default STDOUT, can be changed with output argument
    :param output:
    :param query:
    :param ref:
    :return:
    """
    cmd = blastp_sh.format(ref, qry, qry, eval, threads, output)
    if use_grid == 'T':
        job = bsub(cmd, cpus=threads, name='blastp')
        waitjob(job)
    else:
        sh(cmd)
    return 0


# 0 ref
# 1 qry
# 2 output
blastn_sh = r"""
makeblastdb -in {0} -dbtype nucl
blastn  -num_threads {4} -db {0} -query {1} -evalue {3} -outfmt 6 -out {5}
"""


def blastn(ref=None, qry=None, threads=5, eval=1e-5, use_grid='T', output='-'):
    r"""
    :param query:
    :param ref:
    :return:
    """
    cmd = blastp_sh.format(ref, qry, qry, eval, threads, output)
    if use_grid == 'T':
        job = bsub(cmd, cpus=threads, name='blastn')
        waitjob(job)
    else:
        sh(cmd)
    return 0


def blast2bed(bln=None):
    """
    %s bln > bln.bed
    :param bln:
    :return:
    """
    cmd = r"""
awk '{print $1"\t"$7"\t"$8"\t"$2"\t"$12"\t."}' {0} > {0}.bed
""".format(bln)
    sh(cmd)


def filter_reciprocal_best(bln=None):
    qry_best = {}
    qry_line = {}
    ref_best = {}
    highest = defaultdict(int)
    with open(bln) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            mylist = line.rstrip().split()
            qry = mylist[0]
            ref = mylist[1]
            try:
                bitscore = float(mylist[-1])
            except ValueError:
                logger.error(line)
                continue
            if (bitscore > highest[qry]):
                highest[qry] = bitscore
                qry_best[qry] = ref
                qry_line[qry] = line.rstrip()
            if (bitscore > highest[ref]):
                highest[ref] = bitscore
                ref_best[ref] = qry
    for k in qry_best:
        this_best_ref = qry_best[k]
        if (ref_best[this_best_ref] == k):
            # Reciprocal best
            print(qry_line[k])


def extract_top_n_hits(bln=None, eval=1e-5, top_num=10, output='', threads=4):
    """
    A script function like blastall -v and -b:
    If you used to filter top 5 hits with blastall: blastall -v 5 -b 5
    You can run this script on blast file that has no filter before
    :param top_num: default 10
    :param bln: the blast file (outfmt6 or m8)
    :param eval: default 1e-5
    :return:
    """
    if output == '':
        output = "{0}.top{1}".format(bln, top_num)
    cmd = """
tmp={0}
awk '$11 < 1e-05' {0} > ${{tmp%.raw}}
exit
#debug
awk '$11 < 1e-05' {0} > {0}.filter_eval

sort  --parallel={3} -k1,1 -k12,12gr -k11,11g -k3,3gr {0}.filter_eval > {0}.sorted.qry

# The following command is too slow, deprecated
# Then get the top 5 hits for every query: 
# for next in $(cut -f1 {0}.sorted.qry | sort -u); do grep -w -m {1} "$next" {0}.sorted.qry; done > {0}.sorted.qry.top{1}

sort  --parallel={3} -k2,2 -k12,12gr -k11,11g -k3,3gr {0}.filter_eval > {0}.sorted.ref

# Then get the top 5 hits for every query:
# for next in $(cut -f1 {0}.sorted.ref | sort -u); do grep -w -m {1} "$next" {0}.sorted.ref; done > {0}.sorted.ref.top{1}
""".format(bln, top_num, eval, threads)
    sorted_qry = bln + '.sorted.qry'
    sorted_ref = bln + '.sorted.ref'
    sh(cmd)
    #debug
    return 0
    firstn(sorted_qry, field=1, num_extract=top_num, print_to_screen='F', output=sorted_qry + '.top')
    firstn(sorted_ref, field=2, num_extract=top_num, print_to_screen='F', output=sorted_ref + '.top')
    sh("""cat {0}.sorted.qry.top {0}.sorted.ref.top > {2}""".format(bln, top_num, output))


def extract_reciprocal_best_hits(bln=None, eval=1e-5, top_num=1, output=''):
    """
    A script function like blastall -v and -b:
    If you used to filter top 5 hits with blastall: blastall -v 5 -b 5
    You can run this script on blast file that has no filter before
    :param top_num: default 10
    :param bln: the blast file (outfmt6 or m8)
    :param eval: default 1e-5
    :return:
    """
    if output == '':
        output = "{0}.top{1}".format(bln, top_num)
    cmd = """
awk '$11 < 1e-05' {0} > {0}.filter_eval

sort --parallel=4 -k1,1 -k12,12gr -k11,11g -k3,3gr {0}.filter_eval > {0}.sorted.qry

# The following command is too slow, deprecated
# Then get the top 5 hits for every query: 
# for next in $(cut -f1 {0}.sorted.qry | sort -u); do grep -w -m {1} "$next" {0}.sorted.qry; done > {0}.sorted.qry.top{1}

sort --parallel=4 -k2,2 -k12,12gr -k11,11g -k3,3gr {0}.filter_eval > {0}.sorted.ref

# Then get the top 5 hits for every query:
# for next in $(cut -f1 {0}.sorted.ref | sort -u); do grep -w -m {1} "$next" {0}.sorted.ref; done > {0}.sorted.ref.top{1}
""".format(bln, top_num, eval)
    sorted_qry = bln + '.sorted.qry'
    sorted_ref = bln + '.sorted.ref'
    sh(cmd)
    firstn(sorted_qry, field=1, num_extract=top_num, print_to_screen='F', output=sorted_qry + '.top')
    firstn(sorted_ref, field=2, num_extract=top_num, print_to_screen='F', output=sorted_ref + '.top')
    sh("""cut -f1,2 {0}.sorted.qry.top {0}.sorted.ref.top |sed "s/\t/-/" |sort |uniq -c | grep "2 " |awk '{{print $2}}' > {2}""".format(
        bln, top_num, output))


def firstn(file=None, field=1, num_extract=10, print_to_screen='T', output=''):
    """
    extract first n lines.
    Eg:
        python firstn bln field=1 num_extract=10
        Output(STDOUT):
        first 10 lines of each query (on column1)
    If you want to extract first 10 lines of ref genes of a blast (if it's sorted, it will be top ten hits of ref gene)
        python firstn bln field=2 num_extract=10
    :param file:
    :param field:
    :param num_extract:
    :param print_to_screen: [T/F]
    :return:
    """
    field = int(field)
    num_extract = int(num_extract)
    field = int(field) - 1
    count_dict = defaultdict(int)
    buffer = ''
    with open(file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            count_dict[mylist[field]] += 1
            if count_dict[mylist[field]] > num_extract:
                continue
            else:
                if (print_to_screen == 'T'):
                    print(line, end='')
                buffer += line
    if (output != ''):
        with open(output, 'w') as fh:
            fh.write(buffer)
    else:
        return buffer
    return 0


def filter_bln(bln=None, eval=1e-5, bitscore=0):
    """
    filter blast based on evalue (1e-5 by default) and bitscore (no filter by default)
    output: STDOUT
    :param bln: the fmt6 blast file
    :param eval: float format
    :param bitscore: int format
    :return:
    """
    # cmd = """awk '$11 < {0} && $12 > {1}' {2} """.format(eval, bitscore, bln)
    # result = sh(cmd)
    # print(result)
    eval = float(eval)
    bitscore = int(bitscore)

    with open(bln) as fh:
        for line in fh:
            bitscore_true = False
            eval_true = False

            if line.startswith("#"):
                continue
            mylist = line.rstrip().split()
            try:
                this_bitscore = float(mylist[-1])
                this_eval = float(mylist[-2])
            except ValueError:
                logger.error(line)
                continue
            if this_bitscore > bitscore:
                bitscore_true = True
            if this_eval < eval:
                eval_true = True
            if bitscore_true and eval_true:
                print(line, end='')


# 2021-11-15
# Modified from https://hub.fastgit.org/biopython/biopython/blob/master/Bio/Blast/ParseBlastTable.py
# Fields: query id (1), subject id(2), % identity(3), alignment length(4),
#    mismatches(5), gap opens(6), q. start(7), q. end(8),
#    s. start(9), s. end(10), evalue(11), bit score(12)
class BlastTable:
    """
    Container for Blast Table Entry, the field values from the table.
    qry_range is a list containing two element: query start and query end

    """

    def __init__(self, in_rec):
        """Initialize the class."""
        bt_fields = in_rec.rstrip().split()
        self.qry_id = bt_fields[0]
        self.ref_id = bt_fields[1]
        self.perc_identity = float(bt_fields[2])
        self.aln_len = int(bt_fields[3])
        self.mismatch = int(bt_fields[4])
        self.gaps = int(bt_fields[5])
        self.qry_range = (int(bt_fields[6]), int(bt_fields[7]))
        self.ref_range = (int(bt_fields[8]), int(bt_fields[9]))
        self.e_value = float(bt_fields[10])
        self.bit_score = float(bt_fields[11])

    def get_line(self):
        """
        :return: the formated string (m8 format/outfmt6) of blast object, no return
        """
        mylist = [self.qry_id, self.ref_id,
                  str(self.perc_identity),
                  str(self.aln_len),
                  str(self.mismatch),
                  str(self.gaps),
                  str(self.qry_range[0]), str(self.qry_range[1]),
                  str(self.ref_range[0]), str(self.ref_range[1]),
                  str(self.e_value), str(self.bit_score)]
        return "\t".join(mylist)


if __name__ == "__main__":
    emain()
