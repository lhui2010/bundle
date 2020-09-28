#!/bin/bash
#BSUB -J mergetrain      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail 

##Time Benchmark: 40min

MYSPECIES=kendao
#species=rice
#augustus_model=${species}v1.1
#latin="Oryza sativa"
ROOT=$PWD
REF=$ROOT/input/ref.fa

VERSION=4.0
DIR=$ROOT/step2_annotation/round$VERSION

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

#################
#Train de novo predictor Augustus
#a. Get fasta of gene for training
mkdir -p $DIR/train_augustus 
cd $DIR/train_augustus
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' genome.all.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi $REF -bed - -fo total_v${VERSION}.all.maker.transcripts1000.fasta

cp  /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/../busco/myconfig.ini  ./config.ini
export BUSCO_CONFIG_FILE=$PWD/config.ini
export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

MYSPECIES=kendao
LINEAGE=viridiplantae_odb10
THREADS=104
INPUT=total_v${VERSION}.all.maker.transcripts1000.fasta
OUTPUT=${MYSPECIES}_v${VERSION}
NEWMODEL=${MYSPECIES}_v${VERSION}
AUGUSTUS_SPECIES=rice
AUGUSTUS_CONFIG_PATH_ORIGINAL=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config

if [ -d $OUTPUT ] 
then
    rm -rf $OUTPUT
fi
busco -i $INPUT  -o $OUTPUT  -l $LINEAGE \
  -m genome -c $THREADS --long --augustus_species $AUGUSTUS_SPECIES --augustus_parameters='--progress=true' >busco.out 2>busco.err

cd $OUTPUT/run_viridiplantae_odb10/augustus_output/retraining_parameters/BUSCO_${OUTPUT}/

#kendao_v1.1_maker/run_viridiplantae_odb10/augustus_output/retraining_parameters/BUSCO_kendao_v1.1_maker/
#BUSCO_kendao_v1.1_maker_metapars.cfg
#BUSCO_kendao_v1.1_maker_metapars.cgp.cfg
#BUSCO_kendao_v1.1_maker_metapars.utr.cfg
#BUSCO_kendao_v1.1_maker_parameters.cfg
#BUSCO_kendao_v1.1_maker_weightmatrix.txt

#rename -v 'BUSCO_kendao_v${VERSION}_maker_2982914939' 'kendao_v${VERSION}' *
#perl -e 'for my $f(@ARGV){$new_name=$f; $new_name=~s/BUSCO_kendao_v${VERSION}_maker_\d+/kendao_v${VERSION}/; `mv $f $new_name`;} ' *

rename "BUSCO_" "" *

#sed -i 's/BUSCO_kendao_v${VERSION}_maker_2982914939/kendao_v${VERSION}/g' kendao_v${VERSION}_parameters.cfg
#sed -i 's/BUSCO_kendao_v${VERSION}_maker_[0-9]\+/kendao_v${VERSION}/g' kendao_v${VERSION}_parameters.cfg
sed -i 's/BUSCO_//g' ${MYSPECIES}_v${VERSION}_parameters.cfg

if [ -d $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL ]
then
    RND=$(date +%s%N)
    rm -r $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL 
    #mv $AUGUSTUS_CONFIG_PATH/species/$NEWMODEL $AUGUSTUS_CONFIG_PATH/species/${NEWMODEL}.$RND
fi
mkdir -p $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/${NEWMODEL}
cp ./${OUTPUT}*  $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/${MYSPECIES}_v${VERSION}/
#
echo "Train Augustus completed succefully"
echo "Augustus species: ${MYSPECIES}_v${VERSION}"
date
#################


