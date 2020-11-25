"""
maker relevant utils
"""
import argparse
import re
import sys
from collections import defaultdict
from parse import parse
#from optparse import OptionParser

from iga.apps.base import ActionDispatcher, sh, conda_act, workdir_sh, logger

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
mkdir ${{WORKDIR}}
ccs {0} ${{WORKDIR}}/ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.9
cd ${{WORKDIR}}
samtools view ccs.bam | awk '{{print ">"$1"\n"$10}}' > ccs.fa
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
    isoseq3 polish unpolished.${{i}}.bam ${{ROOT}}/input/*.subreads.bam polished.${{i}}.bam --verbose &
done
wait
samtools merge -@ 20 polished_total.bam polished.*.bam
isoseq3 summarize polished_total.bam summary.csv

# Result is polished_total.bam.fastq 

"""


def isoseq(subreads=None, workdir=''):
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


# 0 ref genome
# 1 qry fasta
fastq2gff_sh = """
minimap2 -t20 -C5 -ax splice {0} {1} |samtools view -F 256 -b - \
| bedtools bamtobed -split -i - > {1}.bed 
gt bed_to_gff3 {1}.bed | sort -k9,9 -k1,1 -k7,7 -k4,4n > {1}.rawgff
python add_match.py ${i}.rawgff > {1}.gff 2> {1}.intron_error
"""


def fastq2gff(fastq=None, genome=None, workdir=''):
    r"""
    align fastqs generated by Trinity or Isoseq to
    references and transform to maker acceptable gffs
    """
    if (workdir == ''):
        workdir = "workdir_isoseq_" + fastq
    cmd = conda_act.format('EDTA') + workdir_sh.format(workdir) + \
          fastq2gff_sh.format(genome, fastq)
    sh(cmd)
    rawgff = fastq + '.rawgff'
    gff = format_gt_gff_to_maker_gff(rawgff)
    return gff


def format_gt_gff_to_maker_gff(gff=None):
    """
    Format gff generated by gt to maker acceptable format
    Input Example:
    000034F|arrow_np1212    .       BED_feature     1899770 1900081 60      -       .       Name=TRINITY_GG_10000_c0_g1_i1
    000034F|arrow_np1212    .       BED_feature     2360009 2360227 60      -       .       Name=TRINITY_GG_10001_c0_g1_i1
    000034F|arrow_np1212    .       BED_feature     2359766 2359901 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    000034F|arrow_np1212    .       BED_feature     2360103 2360189 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    000034F|arrow_np1212    .       BED_feature     2360284 2360402 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    Output Example:
    000034F|arrow_np1212    match   Trinity_Minimap 1899770 1900081 0       -       .       ID=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0
    000034F|arrow_np1212    match_part      Trinity_Minimap 1899770 1900081 60      -       .       ID=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0-exon-1;Name=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0;Parent=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0
    000034F|arrow_np1212    match   Trinity_Minimap 2360009 2360227 0       -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0
    000034F|arrow_np1212    match_part      Trinity_Minimap 2360009 2360227 60      -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0-exon-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0;Parent=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0
    000034F|arrow_np1212    match   Trinity_Minimap 2359766 2361129 0       -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g2_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g2_i1-+-0
    """
    qry1_file = gff
    output_file = gff + ".gff"
    output_buff = ""
    prefix = re.sub(r'\..*', '', qry1_file)
    intron_cutoff = 20000
    source = 'Trinity_Minimap'
    last_chr = ''
    last_end = ''
    last_strand = ''
    last_pos = []
    last_lines = []
    last_name = ''
    count = 0
    feat_count = defaultdict(int)
    with open(qry1_file) as fh:
        for line in fh:
            if (line.startswith('#')):
                continue
            mylist = line.rstrip().split()
            # print(mylist[-1])
            mylist[-1] = mylist[-1].replace('/', '_')
            # Name feats
            this_chr = mylist[0]
            this_type = 'match_part'
            this_start = mylist[3]
            this_end = mylist[4]
            this_score = mylist[5]
            this_strand = mylist[6]
            this_phase = mylist[7]
            feat_parse = parse("Name={name}", mylist[-1])
            try:
                feat_name = prefix + "-+-" + feat_parse['name'] + "-+-" + str(
                    feat_count[feat_parse['name']])  # + this_chr.replace('|', '') + this_start
            except TypeError:
                logger.warning("TypeError on feat {}".format(mylist[-1]))
            this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
            line = "\t".join(
                [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])

            if (last_chr == '' or
                    last_name == feat_name and
                    last_chr == this_chr and
                    last_strand == this_strand and
                    abs(last_end - int(this_start)) < intron_cutoff):
                count += 1

                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join(
                    [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_chr = this_chr
                last_end = int(this_end)
                last_pos.append(int(this_start))
                last_pos.append(int(this_end))
                last_strand = this_strand
                last_lines.append(line)
                last_name = feat_name

            else:
                if (last_name == feat_name and last_chr == this_chr and last_strand == this_strand):
                    logger.warning("LargeIntron {} on {}".format(str(abs(last_end - int(this_start))), feat_name))
                # Prepare print
                match_chr = last_chr
                match_type = "match"
                match_start = str(min(last_pos))
                match_end = str(max(last_pos))
                match_score = '0'
                match_strand = last_strand
                match_phase = '.'
                match_feat = "ID={};Name={}".format(last_name, last_name)
                match_line = "\t".join(
                    [match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase,
                     match_feat])
                last_lines.insert(0, match_line)
                #print("\n".join(last_lines))
                output_buff += "\n".join(last_lines) + "\n"
                # Initialize
                count = 0
                last_chr = this_chr
                last_end = int(this_end)
                last_pos = [int(this_start), int(this_end)]
                last_strand = this_strand
                # Add count of same name so no duplicate name would occur
                feat_name_slim = last_name.split('-+-')[1]
                # exit()
                feat_count[feat_name_slim] += 1
                # print(feat_name_slim)
                # print(feat_count[feat_name_slim])
                try:
                    feat_name = prefix + "-+-" + feat_parse['name'] + "-+-" + str(feat_count[feat_parse['name']])
                except TypeError:
                    logger.warning("TypeError on feat {}".format(mylist[-1]))
                # print(feat_parse['name'])
                # print(feat_count[feat_parse['name']])
                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join(
                    [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_lines = [line]
                last_name = feat_name
        else:
            match_chr = last_chr
            match_type = "match"
            match_start = str(min(last_pos))
            match_end = str(max(last_pos))
            match_score = '0'
            match_strand = last_strand
            match_phase = '.'
            match_feat = "ID={};Name={}".format(last_name, last_name)
            match_line = "\t".join(
                [match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase,
                 match_feat])
            last_lines.insert(0, match_line)
            #print("\n".join(last_lines))
            output_buff += "\n".join(last_lines) + "\n"
    with open(output_file, 'w') as fh:
        fh.write(output_buff)
    return output_file




def str_to_class(str):
    return getattr(sys.modules[__name__], str)


#TODO: use it in the main
def emain(func_name, args):
    # parser = argparse.ArgumentParser(
    #     prog=prog_name,
    #     formatter_class=argparse.RawDescriptionHelpFormatter,
    #     description=textwrap.dedent(usage),
    #     epilog="")
    # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    # args = parser.parse_args()

    # p = argparse.ArgumentParser(prog=func_name, usage=func_doc)
    position_arg = []
    keyword_arg = {}
    number_args = 1
    import sys
    import inspect
    # 下面两行命令用于在函数内部得到函数的名称和文档
    # func_name = sys._getframe().f_code.co_name
    # func_doc = sys._getframe().f_code.co_consts[0]
    # 下面命令用于把字符串的函数名称转换成对象
    #func_name = 'isoseq_'
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

    real_arg = p.parse_args(args)

    # Results for storing arguments after running parse_args
    position_result = []
    keyword_result = {}

    for k in keyword_result:
        keyword_result.update(getattr(real_arg, k))

    for k in position_arg:
        position_result.append(getattr(real_arg, k))

    object_pointer(*position_result, **keyword_result)

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
    # actions = (
    #     ('isoseq', 'extract isoseq flnc reads from subreads.bam')
    #     ('fastq2gff', 'map fastq to reference genome and get gff files'),
    # )
    actions = ['isoseq', 'fastq2gff']
    if(sys.argv[1] in actions):
        action = sys.argv[1]
        if(len(sys.argv) > 2):
            args = sys.argv[2:]
        else:
            args = []
        emain(action, args)
    else:
        print('Possible actions:{}'.format('\n'.join(actions)))
    #p = ActionDispatcher(actions)
    #p.dispatch(globals())

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
