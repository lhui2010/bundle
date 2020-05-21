#!/bin/bash
#
#$ -N train_soft
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -l h=fat01
#

##Time Benchmark: 40min

species=maize
#augustus_model=maize5
latin="Zea mays"
ROOT=$PWD
REF=$ROOT/input/ref.fa

ROUND1=$ROOT/step2_annotation/round1
EST_GFF="$ROUND1/extract_gff/genome.all.noseq.gff.est2genome.gff"
PEP_GFF="$ROUND1/extract_gff/genome.all.noseq.gff.protein2genome.gff"

ROUND_PREVIOUS=$ROOT/step2_annotation/round2

SNAP_HMM=${ROUND_PREVIOUS}/train_snap/snap_round2.hmm
RM_GFF=${ROOT}/step1_repeatmask/Full_mask/full_mask.complex.reformat.gff3
AUGUSTUS_SPECIES=A188_round1

OPT=round3_maker_opts.ctl
##generating CTL files
cd $ROOT/step2_annotation

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
model_org=viridiplantae #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/lustre/local/packages/Maker/maker3/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
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
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
#TMP= #specify a directory other than the system default temporary directory for temporary files
" >${OPT}

#RND=$(date +%s%N)
mkdir round3
cd round3

#25 alway cause error when repeat mask turns on, so we will run batch jobs without 25
for i in `cat ../fa_list.pa`
    do
        echo "#!/bin/bash
#
#$ -cwd
#$ -V
#$ -S /bin/bash
#
set -euxo pipefail

cd $ROOT/step2_annotation/round3
mkdir -p /tmp/\$HOSTNAME/maker.$i
ln -s /tmp/\$HOSTNAME/maker.$i
cd maker.$i
ln -s $ROOT/step2_annotation/fasta_partition/$i.fa
cp $ROOT/step2_annotation/${OPT} ./
cp $ROOT/step2_annotation/maker_exe.ctl ./
cp $ROOT/step2_annotation/maker_bopts.ctl ./
/lustre/local/packages/Maker/maker-3.01.02/bin/maker -g $i.fa $OPT maker_bopts.ctl maker_exe.ctl
cd $ROOT/step2_annotation/round3
mv maker.$i maker_local.$i
cp -r /tmp/\$HOSTNAME/maker.$i maker.$i
" > maker.$i.sh
        qsub maker.$i.sh
        sleep 20s
    done

for i in 25
    do
        echo "#!/bin/bash
#
#$ -cwd
#$ -V
#$ -S /bin/bash
#
set -euxo pipefail

cd $ROOT/step2_annotation/round3
mkdir -p /tmp/\$HOSTNAME/maker.$i
ln -s /tmp/\$HOSTNAME/maker.$i
cd maker.$i
ln -s $ROOT/step2_annotation/fasta_partition/$i.fa
cp $ROOT/step2_annotation/${OPT} ./
cp $ROOT/step2_annotation/maker_exe.ctl ./
cp $ROOT/step2_annotation/maker_bopts.ctl ./
/lustre/local/packages/Maker/maker-3.01.02/bin/maker -R -g $i.fa $OPT maker_bopts.ctl maker_exe.ctl
cd $ROOT/step2_annotation/round3
mv maker.$i maker_local.$i
cp -r /tmp/\$HOSTNAME/maker.$i maker.$i
" > maker.$i.sh
        qsub maker.$i.sh
        sleep 20s
    done
