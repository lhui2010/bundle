#!/usr/bin/env python

import argparse
import textwrap
import subprocess

# class CTL():

# Give arguments, will return the specific CTLs
# Support Direct Print
# Support tag value change
from iga.apps.base import conda_act

nextdenovo_sh = r"""
export PATH=/ds3200_1/users_root/yitingshuang/lh/anaconda2/bin:$PATH

export PATH=$PWD/NextPolish/bin:$PATH
export PATH=$PWD/NextPolish/:$PATH

    WORKDIR="workdir_"${PREFIX}

#ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn

    export PATH=/ds3200_1/users_root/yitingshuang/lh/anaconda2/bin:$PATH

    mkdir -p ${WORKDIR}

#cp NextDenovo/doc/run.cfg ${WORKDIR}/

    cd ${WORKDIR}

 #./01_rundir/03.ctg_graph/01.ctg_graph.sh.work/ctg_graph00/nextgraph.assembly.contig.fasta
         ln -s /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/falcon/workdir_gcpp_on_falconv1/output_resume.fasta coriaria_falcon.fasta
         PREFIX=sgs_polish
         contig=coriaria_falcon.fasta


         touch sgs.fofn ; rm sgs.fofn
         touch ${PREFIX}.0.fa ; rm ${PREFIX}.0.fa
 #        ln -s ${ROOT}/input/${SAMPLE}.sgs.fofn sgs.fofn
         ls /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/clean_fastq/coriaria_sgs_newadp_1.clean.fq.gz  /ds3200_1/users_root/yitingshuang/lh/projects/
 buzzo/clean_fastq/coriaria_sgs_newadp_2.clean.fq.gz> sgs.fofn
         ln -s ${contig} ${PREFIX}.0.fa


         for i in `seq 0 1`
 #        WORKDIR=./03_polish.sgs_round0
 #        cat ${WORKDIR}/03.kmer_count/05.polish.ref.sh.work/polish_genome*/*.fasta > ${PREFIX}.1.fa
 #        for i in `seq 1 2`
         do
             j=`expr $i + 1`
             GENOME=${PREFIX}.${i}.fa
             OUTPUT=${PREFIX}.${j}.fa
             WORKDIR=./03_polish.sgs_round${i}
            echo "[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 30
multithread_jobs = 20
genome = ${GENOME}
genome_size = auto
workdir = ${WORKDIR}
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa
">polish_sgs.cfg
    #echo -e "task = 1212\ngenome = $genome\nsgs_fofn = sgs.fofn" > run.cfg
#Run
            nextPolish polish_sgs.cfg
           # cat ${WORKDIR}/00.sgs_polish/04.polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > ${OUTPUT}
            cat ${WORKDIR}/03.kmer_count/05.polish.ref.sh.work/polish_genome*/*.fasta > ${OUTPUT}
    done


"""

# 0 contig
# 1 bam
# 2 threads
gcpp_sh = """
set -euxo pipefail
pbmm2 align --sort  {0} {1} {0}.bam
samtools index {0}.bam
nproc={2}
gcpp --algorithm=arrow -x 5 -X 120 -q 0 -j $nproc \
        -r {0} {0}.bam \
        -o {0}.polish.fasta,{0}.polish.fastq,{0}.polish.vcf
rm {0}.bam
"""


def gcpp(ctg_file, bam_file, threads=20):
    # assembly, subreads

    cmd = conda_act.format('falcon') + gcpp_sh.format(ctg_file, bam_file, threads)
    #    cmd = conda_act
    subprocess.run(cmd, shell=True)


def main():
    prog_name = "Polish Pacbio Assembly with Pacbio Reads"
    usage = "Polish Pacbio Assembly with Pacbio Reads"

    parser = argparse.ArgumentParser(
        prog=prog_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage),
        epilog="")
    parser.add_argument("CONTIG", help="Raw assembly")
    parser.add_argument("BAM", help="subreads bam")
    parser.add_argument("-t", "--threads", default=100, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()

    ctg_file = args.CONTIG
    bam_file = args.BAM
    gcpp(ctg_file, bam_file)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
