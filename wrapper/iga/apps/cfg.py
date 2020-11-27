from collections import defaultdict

cfg = defaultdict(str)
seperator = defaultdict(str)

__all__ = ['cfg', 'seperator', 'falcon', 'maker']

__author__ = ['Hui Liu <lhui2010@gmail.com>']

__version__ = '0.1'

seperator['maker'] = '='
seperator['falcon'] = '='
cfg['falcon'] = r"""
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

cfg['maker'] = r"""
#-----Genome (these are always required)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

genome=

est= #set of ESTs or assembled mRNA-seq in fasta format
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
rm_gff=${RM_GFF} #pre-identified repeat elements from an external GFF3 file

est_gff=$EST_GFF #aligned ESTs or mRNA-seq from an external GFF3 file
protein_gff=$PEP_GFF  #aligned protein homology evidence from an external GFF3 file

est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no

snaphmm= #SNAP HMM file
augustus_species= #Augustus gene prediction species model


#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
altest= #EST/cDNA sequence file in fasta format from an alternate organism
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)


#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=simple #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=  #  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
gmhmm= #GeneMark HMM file
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no

snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=1 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
#TMP= #specify a directory other than the system default temporary directory for temporary files

"""

cfg['maker_bopts'] = r"""#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'
use_rapsearch=0 #use rapsearch instead of blastx, 1 = yes, 0 = no

pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
"""

cfg['maker_exe'] = r"""
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/ds3200_1/proc/ncbi-blast-2.8.0+/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/ds3200_1/proc/ncbi-blast-2.8.0+/bin/blastn #location of NCBI+ blastn executable
blastx=/ds3200_1/proc/ncbi-blast-2.8.0+/bin/blastx #location of NCBI+ blastx executable
tblastx=/ds3200_1/proc/ncbi-blast-2.8.0+/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch= #location of prerapsearch executable
rapsearch= #location of rapsearch executable
RepeatMasker=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/bin/../exe/RepeatMasker/RepeatMasker #location o
exonerate=/ds3200_1/proc/exonerate-2.2.0-x86_64/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/bin/../exe/snap/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/bin/augustus #loca
fgenesh= #location of fgenesh executable
evm=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/EVidenceModeler-1.1.1/evidence_modeler.pl #location o
tRNAscan-SE=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/bin/../exe/tRNAscan-SE-1.3.1/tRNAscan-SE #locatio
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
"""
