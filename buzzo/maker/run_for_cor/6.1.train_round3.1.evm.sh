#!/bin/bash
#BSUB -J mergetrain      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail 

##Time Benchmark: 40min

MYSPECIES=coriaria
#species=rice
#augustus_model=${species}v1.1
#latin="Oryza sativa"
ROOT=$PWD
REF=$ROOT/input/ref.fa

VERSION=3.1
DIR=$ROOT/step2_annotation.bak/round$VERSION

echo "Start Time:"
date

#################
##First: Merge all fasta files

#a. Create an all-in-one log file
#debug dir
#DIR=$ROOT/step2_annotation.bak/v${VERSION}
#DIR=$ROOT/step2_annotation/round$VERSION
cd $DIR
for i in `cat ../fa_list`;
do
    cat $DIR/maker.$i/$i.maker.output/${i}_master_datastore_index.log |sed "s/\t/\tmaker.$i\/$i.maker.output\//" >>total_master_datastore_index.log
done

#b. Merge fasta
fasta_merge -d total_master_datastore_index.log

#c. Merge GFF, two-step approach
for i in `cat ../fa_list` 
do  
    cd $DIR/maker.${i}/${i}.maker.output && gff3_merge -d ${i}_master_datastore_index.log
done

cd $DIR
gff3_merge  -o $DIR/genome.all.gff $DIR/maker.*/*.maker.output/*.gff
gff3_merge -n -o $DIR/genome.all.noseq.gff $DIR/maker.*/*.maker.output/*.gff

echo "Merge completed succefully:"
date

#################


#################
#Split_GFF files
mkdir -p $DIR/extract_gff
cd $DIR/extract_gff
ln -s ../genome.all.noseq.gff
split_gff.pl genome.all.noseq.gff
echo "est_gff: $DIR/extract_gff/genome.all.noseq.gff.est2genome.gff"
echo "pep_gff: $DIR/extract_gff/genome.all.noseq.gff.protein2genome.gff"
#################

