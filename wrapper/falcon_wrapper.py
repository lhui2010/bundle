#!/usr/bin/env python

import argparse
import textwrap
import subprocess

#class CTL():


cfg = {}
cfg['falcon_cfg']="""
#### Input
[General]
input_fofn=Eu_fasta.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=streamed-median
target=assembly
skip_checks=False
LA4Falcon_preload=false

#### Data Partitioning
pa_DBsplit_option=-x500 -s400
ovlp_DBsplit_option=-s400

#### Repeat Masking
pa_HPCTANmask_option=
#no-op repmask param set
pa_REPmask_code=0,300;0,300;0,300

####Pre-assembly
# adjust to your genome size
genome_size = 671248000
seed_coverage = 40
length_cutoff = -1
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option= -k18 -e0.80 -l1000 -h256 -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 4 --max-n-read 200
falcon_sense_greedy=False

####Pread overlapping
ovlp_HPCdaligner_option=-v -B128 -M24
ovlp_daligner_option= -k24 -e.92 -l1800 -h1024 -s100

####Final Assembly
length_cutoff_pr=1000
overlap_filtering_setting=--max-diff 100 --max-cov 100 --min-cov 2
fc_ovlp_to_graph_option=

[job.defaults]
job_type=lsf
pwatcher_type=blocking
JOB_QUEUE=Q104C512G_X4
MB=32768
NPROC=2
njobs=100
submit = bsub -K -n 2 -q ${JOB_QUEUE} -J ${JOB_NAME} -o ${JOB_STDOUT} -e ${JOB_STDERR} ${JOB_SCRIPT}

[job.step.da]
NPROC=4
MB=32768
njobs=200
[job.step.la]
NPROC=4
MB=32768
njobs=200
[job.step.cns]
NPROC=4
MB=32768
njobs=200
[job.step.pda]
NPROC=4
MB=32768
njobs=200
[job.step.pla]
NPROC=4
MB=32768
njobs=200
[job.step.asm]
NPROC=24
MB=196608
njobs=1
"""

class CTF():
    """
    Give arguments, will return the specific CTLs
    Support Direct Print
    Support tag value change
    """
    def __init__(self, type="falcon"):
        self.content = cfg[type]

    def change_val(self, args):
            template_name="Falcon", parse_args):
        ar

    def printout(self):
        print self.conent

    @static
    def get_fofn(file_list, fofn_file):

    def prepare_workdir():

#0 falcon input
falcon_sh = """
fc_run {0}
"""

def falcon_run(fasta_file, genome_size, threads=20):
#assembly, subreads
    conda_act = r"""
    export PS1="(base) \[\033]2;\h:\u $PWD\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
    source ~/lh/anaconda3/etc/profile.d/conda.sh
    conda activate falcon
    """
    ctf_file = CTF('falcon')
    ctf_file.change_val('input_fofn={0}; genome_size={1}'.format(fasta_file, genome_size)
    ctf_file.print_out('/tmp/output')
    cmd = conda_act + gcpp_sh.format(ctg_file, bam_file, threads)
#    cmd = conda_act 
    subprocess.run(cmd, shell = True)


def main():
    prog_name = "Polish Pacbio Assembly with Pacbio Reads"
    usage = "Polish Pacbio Assembly with Pacbio Reads"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("fasta", help="Raw assembly")
    args = parser.parse_args()  

    get_fofn(args.fasta, 'fofn_path')
    ctg_file = args.CONTIG
    falcon_run(ctg_file, bam_file)
#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
