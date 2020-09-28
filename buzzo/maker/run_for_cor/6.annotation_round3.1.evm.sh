#!/bin/bash
#
#$ -N train_soft
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -l h=fat01
#

##Time Benchmark: 40min

set -euxo pipefail


MYSPECIES=coriaria
ROOT=$PWD
REF=$ROOT/input/ref.fa

PREV_VERSION=2.0
VERSION=3.1
PROJECTID=round${VERSION}

ROUND_PREV=$ROOT/step2_annotation/round${PREV_VERSION}

#Prepare for annotated gffs
#...
EST_GFF=$ROOT/input/transcripts.gff
PEP_GFF=$ROOT/input/peps.gff
#./step2_annotation/round1/extract_gff/genome.all.noseq.gff.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file
#EST_GFF="$ROUND_PREV/extract_gff/genome.all.noseq.gff.est2genome.gff"
#PEP_GFF="$ROUND_PREV/extract_gff/genome.all.noseq.gff.protein2genome.gff"
#/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/step2_annotation/round1.1/train_snap/snap_v1.1.hmm
#/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM//kendao_v1.1.hmm
#kendao_v1.1
#SNAP_HMM=$ROUND_PREV/train_snap/snap_v${PREV_VERSION}.hmm
ZOE_HMM_DIR=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/
SNAP_HMM=${ZOE_HMM_DIR}/${MYSPECIES}_v${PREV_VERSION}.hmm
RM_GFF=${ROOT}/step1_repeatmask/Full_mask/full_mask.complex.reformat.gff3
AUGUSTUS_SPECIES=${MYSPECIES}_v${PREV_VERSION}

##generating CTL files
cd $ROOT/step2_annotation

OPT=round${VERSION}_maker_opts.ctl

echo "
#-----Genome (these are always required)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic


est= #set of ESTs or assembled mRNA-seq in fasta format
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
rm_gff=$RM_GFF #pre-identified repeat elements from an external GFF3 file

est_gff=$EST_GFF #aligned ESTs or mRNA-seq from an external GFF3 file
protein_gff=$PEP_GFF  #aligned protein homology evidence from an external GFF3 file

est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no

snaphmm=$SNAP_HMM #SNAP HMM file
augustus_species=$AUGUSTUS_SPECIES #Augustus gene prediction species model


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
model_org=rice #select a model organism for RepBase masking in RepeatMasker
rmlib=${ROOT}/step1_repeatmask/ref-families.fa  #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
gmhmm= #GeneMark HMM file
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=1 #run EvidenceModeler, 1 = yes, 0 = no

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
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
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
" >${OPT}

#maker -CTL

if [ -d ${PROJECTID} ]
then
    RND=$(date +%s%N)
    mv ${PROJECTID} ${PROJECTID}.$RND
fi

mkdir -p ${PROJECTID}

for i in `cat fa_list`
    do
        cd $ROOT/step2_annotation/${PROJECTID}
        mkdir -p maker.$i
        cd maker.$i
        #mkdir /tmp/
        #mkdir maker.$i
        #cd maker.$i
        ln -s ../../fasta_partition/$i.fa
        cp ../../${OPT} ./
        cp ../../../maker_exe.ctl ./
        cp ../../../maker_bopts.ctl ./
        cp ../../../maker_evm.ctl ./
        echo "genome=$PWD/$i.fa #genome sequence (fasta file or fasta embeded in GFF3 file)" >>$OPT
       # bsub -J mk$i -o maker.$i.log -e maker.$i.err "maker -R $OPT maker_bopts.ctl maker_exe.ctl >>"
        bsub -q Q104C512G_X4 -J R3mk$i -o output.$i -e error.$i " export AUGUSTUS_CONFIG_PATH=/tmp/lh_config && maker $OPT maker_bopts.ctl maker_exe.ctl maker_evm.ctl >> maker.$i.log 2>>maker.$i.err"
        sleep 10s
 #       echo "/lustre/eocal/packages/Maker/maker-3.01.02/bin/maker $OPT maker_bopts.ctl maker_exe.ctl" >run.sh
    done
