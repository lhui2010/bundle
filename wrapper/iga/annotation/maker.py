"""
maker relevant utils
"""
import argparse
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


def isoseq_(subreads, workdir=''):
    if(type(subreads) == list):
        subreads = " ".join(subreads)
    if (workdir == ''):
        workdir = "workdir_isoseq_" + subreads.split()[0]
    cmd = conda_act.format('isoseq3') + isoseq_sh.format(subreads, workdir)
    sh(cmd)

def isoseq(args):
    """
    %prog isoseq subreads.fasta

    Wrapper for `isoseq`
    """
    # parser = argparse.ArgumentParser(
    #     prog=prog_name,
    #     formatter_class=argparse.RawDescriptionHelpFormatter,
    #     description=textwrap.dedent(usage),
    #     epilog="")
    # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    # args = parser.parse_args()
    import sys
    func_name = sys._getframe().f_code.co_name
    p = argparse.ArgumentParser(prog=func_name, usage=__doc__)
    p.add_argument("subreads", help="subreads bams from Isoseq", nargs='+')
    p.add_argument("-d", "--workdir", help="Name of working directory")
    p.add_argument("-o", "--output", help="Output file name")

    p.parse_args(args)

    subreads_files = p.subreads
    isoseq_(subreads_files)


def minimap_rna(transcript, genome, threads=30, output=''):
    if (output == ''):
        output = transcript + ".sam"
    cmd = minimap_rna_sh.format(threads, genome, transcript, output)
    sh(cmd)


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
