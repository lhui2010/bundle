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

PREVVERSION=2.0
VERSION=3.0
DIR=$ROOT/step2_annotation/round$VERSION

echo "Start Time:"
date

#################
##First: Merge all fasta files

#a. Create an all-in-one log file
#debug dir
#DIR=$ROOT/step2_annotation.bak/v${VERSION}
#DIR=$ROOT/step2_annotation/round$VERSION

#################
#Train de novo predictor Augustus
#a. Get fasta of gene for training
mkdir -p $DIR/train_augustus 
cd $DIR/train_augustus
#../extract_gff/genome.all.noseq.gff.maker.gff
if [ ! -e genome.all.gff ]
then
    ln -s ../extract_gff/genome.all.noseq.gff.maker.gff genome.all.gff
fi

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' genome.all.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi $REF -bed - -fo total_v${VERSION}.all.maker.transcripts1000.fasta

cp  /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/../busco/myconfig.ini  ./config.ini
export BUSCO_CONFIG_FILE=$PWD/config.ini
export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

LINEAGE=viridiplantae_odb10
LINEAGE=embryophyta_odb10
THREADS=104
INPUT=total_v${VERSION}.all.maker.transcripts1000.fasta
OUTPUT=${MYSPECIES}_v${VERSION}
NEWMODEL=${MYSPECIES}_v${VERSION}
AUGUSTUS_SPECIES=${MYSPECIES}_v${PREVVERSION}
AUGUSTUS_CONFIG_PATH_ORIGINAL=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config

if [ -d $OUTPUT ] 
then
    rm -rf $OUTPUT
fi
busco -i $INPUT  -o $OUTPUT  -l $LINEAGE \
  -m genome -c $THREADS --long --augustus_species $AUGUSTUS_SPECIES --augustus_parameters='--progress=true' >busco.out 2>busco.err

cd $OUTPUT/run_${LINEAGE}/augustus_output/retraining_parameters/BUSCO_${OUTPUT}/

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


