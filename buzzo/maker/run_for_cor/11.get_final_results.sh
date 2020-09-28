#!/bin/bash
#
#$ -N Get_Results
#$ -cwd
#$ -V
#$ -S /bin/bash
#


species=maize
SPECIES_NAME=A188
augustus_model=A188_round1
#augustus_model=maize5
latin="Zea mays"
ROOT=$PWD
REF=$ROOT/input/ref.fa
RM_GFF=${ROOT}/step1_repeatmask/Full_mask/full_mask.complex.reformat.gff3
MAKER_GFF=${ROOT}/step2_annotation/round3/extract_gff/genome.all.noseq.gff.maker.gff
GFF=genome.all.noseq.gff.maker.gff
SCRIPTS_DIR=${ROOT}/scripts
PEP=${ROOT}/step2_annotation/round3/total.all.maker.proteins.fasta
CDS=${ROOT}/step2_annotation/round3/total.all.maker.transcripts.fasta

echo "Start Time:"
date

#DIR=$ROOT/step3_filter_annotation/2.filter_aed
DIR=$ROOT/step3_filter_annotation/3.results
mkdir -p $DIR && cd $DIR

#wc -l final_list.txt
#Extract finalGFF
perl ${SCRIPTS_DIR}/filter_GFF.pl ../2.filter_aed/final_list.txt ../2.filter_aed/${GFF} >final_list.GFF
#Get Longest transcripts & GeneID table
python ${SCRIPTS_DIR}/longest_transcript_from_gff.py final_list.GFF > final_list.GFF.longest.table
#Extract longest gff
perl ${SCRIPTS_DIR}/filter_GFF2.pl final_list.GFF.longest.table final_list.GFF > final_list.GFF.longest

#Extract longest PEP
perl ${SCRIPTS_DIR}/select_fastav1.pl final_list.GFF.longest.table $PEP >final_list.pep
#Extract longest CDS
perl ${SCRIPTS_DIR}/select_fastav1.pl final_list.GFF.longest.table $CDS >final_list.cds


###Change IDs###
maker_map_ids --prefix ${SPECIES_NAME}G --justify=5 --suffix='-t' --iterate='1' final_list.GFF.longest  >${SPECIES_NAME}.gff.map
cp final_list.GFF.longest ${SPECIES_NAME}.gff
map_gff_ids ${SPECIES_NAME}.gff.map ${SPECIES_NAME}.gff
cp final_list.pep ${SPECIES_NAME}.pep
map_fasta_ids ${SPECIES_NAME}.gff.map ${SPECIES_NAME}.pep
cp final_list.cds ${SPECIES_NAME}.cds
map_fasta_ids ${SPECIES_NAME}.gff.map ${SPECIES_NAME}.cds





