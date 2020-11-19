
from collections import defaultdict
cfg = defaultdict(str)
seperator = defaultdict(str)

__all__=['cfg', 'seperator', 'falcon', 'maker']

__author__ = ['Hui Liu <lhui2010@gmail.com>']

__version__ = '0.1'

seperator['falcon'] = '='
cfg['falcon']=r"""
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
