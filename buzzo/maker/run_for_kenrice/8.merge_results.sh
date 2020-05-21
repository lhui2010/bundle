#!/bin/bash
#
#$ -N merge_round3
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -l h=fat01
#

#fragScaff_scaffold_135是出错后面手动跑的，所以单独修改了脚本加了进来
##Time Benchmark: 40min

species=maize
augustus_model=A188_round1
#augustus_model=maize5
latin="Zea mays"
ROOT=$PWD
REF=$ROOT/input/ref.fa

echo "Start Time:"
date

#################
##First: Merge all fasta files

#a. Create an all-in-one log file
#debug dir
#DIR=$ROOT/step2_annotation.bak/round2
DIR=$ROOT/step2_annotation/round3
cd $DIR
for i in `cat ../fa_list`;
do
    cat $DIR/maker.$i/$i.maker.output/${i}_master_datastore_index.log |sed "s/\t/\tmaker.$i\/$i.maker.output\//" >>total_master_datastore_index.log
done

for i in 135R
    do
        cat $DIR/maker.22/$i.maker.output/${i}_master_datastore_index.log |sed "s/\t/\tmaker.$i\/$i.maker.output\//" >>total_master_datastore_index.log
done

#b. Merge fasta
fasta_merge -d total_master_datastore_index.log

#c. Merge GFF, two-step approach
for i in `cat ../fa_list` 
do  
    cd $DIR/maker.${i}/${i}.maker.output && gff3_merge -d ${i}_master_datastore_index.log
done

for i in 135R
do
    cd $DIR/maker.22/${i}.maker.output && gff3_merge -d ${i}_master_datastore_index.log
done

cd $DIR
gff3_merge  -o $DIR/genome.all.gff $DIR/maker.*/*.maker.output/*.gff
gff3_merge -n -o $DIR/genome.all.noseq.gff $DIR/maker.*/*.maker.output/*.gff

echo "Merge completed succefully:"
date

#################


#################
#Split_GFF files
mkdir $DIR/extract_gff
cd $DIR/extract_gff
ln -s ../genome.all.noseq.gff
split_gff.pl genome.all.noseq.gff
echo "est_gff: $DIR/extract_gff/genome.all.noseq.gff.est2genome.gff"
echo "pep_gff: $DIR/extract_gff/genome.all.noseq.gff.protein2genome.gff"
#################

