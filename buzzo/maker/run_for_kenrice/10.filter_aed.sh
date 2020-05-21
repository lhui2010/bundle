#!/bin/bash
#
#$ -N filter_repeat
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
RM_GFF=${ROOT}/step1_repeatmask/Full_mask/full_mask.complex.reformat.gff3
MAKER_GFF=${ROOT}/step2_annotation/round3/extract_gff/genome.all.noseq.gff.maker.gff
GFF=genome.all.noseq.gff.maker.gff

echo "Start Time:"
date

DIR=$ROOT/step3_filter_annotation/2.filter_aed
mkdir -p $DIR && cd $DIR

ln -s $MAKER_GFF
grep --color=auto -P "\tmRNA\t"  $GFF |awk '{print $9}' | sed 's/ID=//; s/;Parent=/\t/; s/;Name=.*_AED=/\t/; s/;.*$//' >${GFF}.aed_table

awk '$3<=0.2' ${GFF}.aed_table >${GFF}.aed_table.filter

cut -f2 ${GFF}.aed_table.filter |sort |uniq >${GFF}.aed_table.filter.genelist

echo "Number of AED<=0.2 genes"
wc -l ${GFF}.aed_table.filter.genelist

#;_QI=
awk '$3<=0.5' ${GFF}.aed_table >${GFF}.aed_table.round2
#snap-gene-0.0-mRNA-1;_AED=0.24;_eAED=0.37;_QI=0|0|0.33|1|0|0.33|3|231|247
grep --color=auto -P "\tmRNA\t"  $GFF |awk '{print $9}' | sed 's/ID=//; s/;Parent=/\t/; s/;Name=.*_QI=/\t/; s/;.*$//; s/|/\t/g' >${GFF}.QI_table
awk '$4>0 ' ${GFF}.QI_table >${GFF}.QI_table.ss
awk '$6>0.5' ${GFF}.QI_table >${GFF}.QI_table.qhlr
for i in ${GFF}.aed_table.round2 ${GFF}.QI_table.ss  ${GFF}.QI_table.qhlr
    do
        cut -f2 $i |sort |uniq >$i.genelist
        done
wc -l *.genelist
for i in ${GFF}.QI_table.qhlr.genelist ${GFF}.QI_table.ss.genelist ; do selectIterm.pl ${GFF}.aed_table.round2.genelist $i >$i.selected; done

cat *selected ${GFF}.aed_table.filter.genelist |sort |uniq >final_list.txt
echo "Final gene sets:"
wc -l final_list.txt

