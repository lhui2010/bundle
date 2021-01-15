"""
Isoseq relevant utils
"""

from iga.apps.base import sh, emain, conda_act

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

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
isoseq3 cluster isoseq_flnc.bam unpolished.bam --verbose --use-qvs
#isoseq3 cluster isoseq_flnc.bam unpolished.bam --split-bam 10
# pbindex subreads.bam
# for i in `seq 0 9`
# do
#     isoseq3 polish unpolished.${{i}}.bam ${{ROOT}}/input/*.subreads.bam polished.${{i}}.bam --verbose &
# done
# wait
# samtools merge -@ 20 polished_total.bam polished.*.bam
# isoseq3 summarize polished_total.bam summary.csv
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


# 0 workdir
# 1 subreads.bam
# 2 primer.fa
isoseq_pb_sh = r"""mkdir -p {0}
cd {0}
ln -s ../{1}
ln -s ../{2}
ccs {1} {1}.ccs.bam --min-rq 0.9
#lima is used to remove primer sequence
#but can it be used to identify reads containing primer sequence as full length reads?
lima {1}.ccs.bam {2} {1}.fl.bam --isoseq  --peek-guess
#flnc equals full-length, non-concatemer
INPUTBAM={1}
PRIMER={2}
isoseq3 refine ${{INPUTBAM%.bam}}.fl.*.bam ${{PRIMER}} ${{INPUTBAM%.bam}}.flnc.bam --require-polya
isoseq3 cluster ${{INPUTBAM%.bam}}.flnc.bam ${{INPUTBAM%.bam}}.clustered.bam --verbose --use-qvs
"""


def isoseq_pb(subreads=None, workdir=''):
    r"""
    convert Isoseq(pacbio standard) subreads.bam to flnc.fastq
    :param subreads:
    :param workdir:
    :return:
    """
    if (type(subreads) == list):
        subreads = " ".join(subreads)
    if (workdir == ''):
        workdir = "workdir_isoseq_" + subreads.split()[0]
    cmd = conda_act.format('isoseq3') + isoseq_sh.format(subreads, workdir)
    sh(cmd)

if __name__ == "__main__":
    emain()