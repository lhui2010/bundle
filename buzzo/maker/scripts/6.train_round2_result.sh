#!/bin/bash
#
#$ -N train_soft_round2
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -l h=fat01
#

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
DIR=$ROOT/step2_annotation/round2
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
mkdir $DIR/extract_gff
cd $DIR/extract_gff
ln -s ../genome.all.noseq.gff
split_gff.pl genome.all.noseq.gff
echo "est_gff: $DIR/extract_gff/genome.all.noseq.gff.est2genome.gff"
echo "pep_gff: $DIR/extract_gff/genome.all.noseq.gff.protein2genome.gff"
#################

#################
#Train de novo predictor SNAP
mkdir $DIR/train_snap
cd $DIR/train_snap
ln -s ../genome.all.gff
maker2zff -x 0.25 -l 50  genome.all.gff
fathom -gene-stats genome.ann genome.dna >gene-stats.log 2>&1
fathom -validate genome.ann genome.dna >validate.log 2>&1

fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna

mkdir params

cd params

forge ../export.ann ../export.dna >../forge.log 2>&1

cd ..
hmm-assembler.pl snap_round2 params > snap_round2.hmm

echo "Train SNAP completed succefully:"
echo "$DIR/train_snap/snap_round2.hmm"
date
#################
#Train de novo predictor
#a. Get fasta of gene for training
mkdir $DIR/train_augustus 
cd $DIR/train_augustus
ln -s ../genome.all.gff

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' genome.all.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi $REF -bed - -fo total_rnd1.all.maker.transcripts1000.fasta

export PYTHONPATH=/lustre/local/lib/python2.7/site-packages/:\$PYTHONPATH
cp /lustre/home/liuhui/bin/busco_config.ini ./config.ini
export BUSCO_CONFIG_FILE=$PWD/config.ini

run_BUSCO.py -i total_rnd1.all.maker.transcripts1000.fasta  -o A188_rnd1_maker -l /lustre/local/database/BUSCODB/embryophyta_odb9/ \
  -m genome -c 160 --long -sp $augustus_model -z --augustus_parameters='--progress=true' >busco.out 2>busco.err

cd run_A188_rnd1_maker/augustus_output/retraining_parameters

#rename -v 'BUSCO_A188_rnd1_maker_2982914939' 'A188_round2' *
perl -e 'for my $f(@ARGV){$new_name=$f; $new_name=~s/BUSCO_A188_rnd1_maker_\d+/A188_round2/; `mv $f $new_name`;} ' *

#sed -i 's/BUSCO_A188_rnd1_maker_2982914939/A188_round2/g' A188_round2_parameters.cfg
sed -i 's/BUSCO_A188_rnd1_maker_[0-9]\+/A188_round2/g' A188_round2_parameters.cfg

RND=$(date +%s%N)
mv $AUGUSTUS_CONFIG_PATH/species/A188_round2 $AUGUSTUS_CONFIG_PATH/species/A188_round2.$RND
mkdir $AUGUSTUS_CONFIG_PATH/species/A188_round2
cp A188_round2*  $AUGUSTUS_CONFIG_PATH/species/A188_round2/
#
echo "Train Augustus completed succefully"
echo "Augustus species: A188_round2"
date
#################



