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

echo "Start Time:"
date

DIR=$ROOT/step3_filter_annotation/1.filter_repeats
mkdir -p $DIR && cd $DIR

ln -s $RM_GFF
ln -s $MAKER_GFF
TOTALGFF=`basename $MAKER_GFF`
GENEGFF=${TOTALGFF%.*}.gene.gff
grep -P "\tgene\t" $TOTALGFF >$GENEGFF

#bedtools intersect -b $RM_GFF -a $GENEGFF -wo |perl -e 'while(<>){chomp; my @e=split/\t/, $_; $per=$e[-1]/($e[4]-$e[3]); print $_."\t$per\n";}' >maker_gene_repeat_overlap.gff

bedtools coverage -a $GENEGFF -b $RM_GFF  >${GENEGFF}.overlap.gff

#cut -f20 maker_gene_repeat_overlap.gff >maker_gene_repeat_overlap.gff.dist

echo "a=read.table(\"${GENEGFF}.overlap.gff\")
pdf(\"dist.pdf\")
hist(a[,13])
" >hist.R
Rscript hist.R
#awk -F"\t" '$20<0.5' maker_gene_repeat_overlap.gff >maker_gene_repeat_overlap.gff.filter

awk '$13<0.5' ${GENEGFF}.overlap.gff >${GENEGFF}.overlap.filter.gff

echo "Number of Filtered genes that overlap with TE <50% "
wc -l ${GENEGFF}.overlap.filter.gff


