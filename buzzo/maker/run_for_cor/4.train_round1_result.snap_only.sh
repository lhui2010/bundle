#!/bin/bash
#BSUB -J merge_snap      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail 

##Time Benchmark: 40min

ROOT=$PWD
REF=$ROOT/input/ref.fa

VERSION=1

echo "Start Time:"
date

DIR=$ROOT/step2_annotation/round$VERSION
cd $DIR

MYSPECIES=coriaria
THREADS=104
INPUT=total_v${VERSION}.all.maker.transcripts1000.fasta
OUTPUT=${MYSPECIES}_v${VERSION}
NEWMODEL=${MYSPECIES}_v${VERSION}


echo "Start Time:"
date

#################
#Train de novo predictor SNAP
HMMDIR=~/lh/bin/maker3/exe/snap/Zoe/HMM/
mkdir -p $DIR/train_snap
cd $DIR/train_snap
ln -s ../genome.all.gff
maker2zff -x 0.25 -l 50  genome.all.gff
fathom -gene-stats genome.ann genome.dna >gene-stats.log 2>&1
fathom -validate genome.ann genome.dna >validate.log 2>&1

fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna

mkdir -p params

cd params

forge ../export.ann ../export.dna >../forge.log 2>&1

cd ..
hmm-assembler.pl snap_v${VERSION} params > snap_v${VERSION}.hmm

cp snap_v${VERSION}.hmm ${HMMDIR}/${MYSPECIES}_v${VERSION}.hmm

echo "Train SNAP completed succefully:"
echo "$DIR/train_snap/snap_v${VERSION}.hmm"
echo "${HMMDIR}/${MYSPECIES}_v${VERSION}.hmm"
date
###################
