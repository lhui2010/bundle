"""
maker relevant utils
"""
import argparse
import sys
from optparse import OptionParser

from iga.apps.base import ActionDispatcher, sh, conda_act

# def sam2gff(sam, gff=""):
#
# isoseq_official_sh=r"""
# ccs [movie].subreads.bam [movie].ccs.bam --min-rq 0.9
# lima --isoseq --dump-clips --no-pbi --peek-guess -j 24 ccs.bam primers.fasta demux.bam
# isoseq3 refine --require-polya combined_demux.consensusreadset.xml primers.fasta flnc.bam
# bamtools convert -format fastq -in flnc.bam > flnc.fastq
# """

# 0 subreads.bam
# 1 workdir
# 2 output.bam
isoseq_sh = r"""export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/isoseq3/BGI-Full-Length-RNA-Analysis-Pipeline/bin:$PATH
export PERL5LIB=""
WORKDIR={1}
ccs {0} ${WORKDIR}/ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.75
cd ${WORKDIR}
samtools view ccs.bam | awk '{print ">"$1"\n"$10}' > ccs.fa
echo ">primer_F
AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG
>primer_S
GTACTCTGCGTTGATACCACTGCTTACTAGT">primer.fa
makeblastdb -in primer.fa -dbtype nucl
blastn -num_threads 32 -query ccs.fa -db primer.fa -outfmt 7 -word_size 5 > mapped.m7
classify_by_primer -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 16 -min_isolen 200 -outdir ./
flnc2sam ccs.sam isoseq_flnc.fasta > isoseq_flnc.sam
samtools view -bS isoseq_flnc.sam > isoseq_flnc.bam
isoseq3 cluster isoseq_flnc.bam unpolished.bam --split-bam 10
pbindex subreads.bam
for i in `seq 0 9`
do
    isoseq3 polish unpolished.${i}.bam ${ROOT}/input/*.subreads.bam polished.${i}.bam --verbose &
done
wait
samtools merge -@ 20 polished_total.bam polished.*.bam
isoseq3 summarize polished_total.bam summary.csv

# Result is polished_total.bam.fastq 

"""


def isoseq_(subreads=None, workdir=''):
    r"""
    isoseq subreads.fasta

    Wrapper for `isoseq`
    """
    if (type(subreads) == list):
        subreads = " ".join(subreads)
    if (workdir == ''):
        workdir = "workdir_isoseq_" + subreads.split()[0]
    cmd = conda_act.format('isoseq3') + isoseq_sh.format(subreads, workdir)
    sh(cmd)


def str_to_class(str):
    return getattr(sys.modules[__name__], str)


def isoseq(args):

    # parser = argparse.ArgumentParser(
    #     prog=prog_name,
    #     formatter_class=argparse.RawDescriptionHelpFormatter,
    #     description=textwrap.dedent(usage),
    #     epilog="")
    # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    # args = parser.parse_args()

    #p = argparse.ArgumentParser(prog=func_name, usage=func_doc)
    position_arg = []
    keyword_arg = {}
    number_args = 1
    import sys
    import inspect
    # 下面两行命令用于在函数内部得到函数的名称和文档
    # func_name = sys._getframe().f_code.co_name
    # func_doc = sys._getframe().f_code.co_consts[0]
    # 下面命令用于把字符串的函数名称转换成对象
    func_name = 'isoseq_'
    object_pointer = getattr(sys.modules[__name__], func_name)
    p = argparse.ArgumentParser(prog=func_name, usage=object_pointer.__doc__)
    # 下面的两个命令用于从函数对象中调取形参的名字和默认值（空值用Nonetype表示），用来转换成parse_args
    for kw, kw_defaults in zip(inspect.getfullargspec(object_pointer).args,
                               inspect.getfullargspec(object_pointer).defaults):
        if (kw_defaults == None):
            position_arg.append(kw)
        else:
            keyword_arg[kw] = kw_defaults
    if (len(position_arg) == 1):
        # If only one input arg is needed for the function, allow multiple files as input
        number_args = '+'
    for k in position_arg:
        p.add_argument(k, help=k, nargs=number_args)
    for k, v in keyword_arg.items():
        p.add_argument("--" + k, default=v)

    p.parse_args(args)

    print(dir(p))
    print(p.subreads)

#Results for storing arguments after running parse_args
    position_result = []
    keyword_result = {}

    for k in keyword_result:
        keyword_result.update(getattr(p, k))

    for k in position_arg:
        position_result.append(getattr(p, k))

    object_pointer(**position_result, **keyword_result)

    # if(number_args == 1):
    #     position_result = getattr(p, position_arg[0])
    #     object_pointer(**position_result, **keyword_result)
    # else:
    #     for k in position_arg:
    #         position_result.append(getattr(p, k))
    # object_pointer(p.__dict__)

# def minimap_rna(transcript, genome, threads=30, output=''):
#     if (output == ''):
#         output = transcript + ".sam"
#     cmd = minimap_rna_sh.format(threads, genome, transcript, output)
#     sh(cmd)


# parallel run

# train

# liftover
# Require RaGOO


def liftover():
    pass


# function

def main():
    """
    """
    actions = (
        ('isoseq', 'extract isoseq flnc reads from subreads.bam'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

    # print(__file__)
    # print(__doc__)
    # exit()
    # prog_name = "busco_wrapper"
    # usage = "run busco on selected GENOME"
    #
    # parser = argparse.ArgumentParser(
    #     prog=prog_name,
    #     formatter_class=argparse.RawDescriptionHelpFormatter,
    #     description=textwrap.dedent(usage),
    #     epilog="")
    # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    # args = parser.parse_args()
    #
    # busco(args.GENOME)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
