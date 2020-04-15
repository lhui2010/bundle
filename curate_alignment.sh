#!/bin/bash

GENE_ID=$1
#M445T00004237

if [ -f ${GENE_ID}.combine.fa.aln ];then
    echo "Already  exist ${GENE_ID}, exiting.."
    exit
fi
CHR=`grep ${GENE_ID} annotation/M445.bed|awk '{print $1}'`
START=`grep ${GENE_ID} annotation/M445.bed|awk '{print $2}'`
START=$((START - 2000))
END=`grep ${GENE_ID} annotation/M445.bed|awk '{print $3}'`
END=$((END + 2000))
QRY_LEN=`expr ${END} - ${START}`

#echo $CHR
#echo $START
#echo $END
#
#exit
fa_pos.pl M445-4.chr.fa ${CHR} ${START} ${END}  > ${GENE_ID}.gene.fa
lastal -u 0 -P 56 -i3G -f BlastTab M441-5.chr ${GENE_ID}.gene.fa | grep -v "^#"> ${GENE_ID}.gene.fa.last
ORTHO_CHR=`head -1 ${GENE_ID}.gene.fa.last | awk '{print $2}'`
ORTHO_START=`head -1 ${GENE_ID}.gene.fa.last | awk '{print $9}'`
ORTHO_END=`head -1 ${GENE_ID}.gene.fa.last | awk '{print $10}'`

ORTHO_LEN=`expr ${ORTHO_END} - ${ORTHO_START}`

if [ ${QRY_LEN} -ne ${ORTHO_LEN} ]; then
    echo "Error when processing gene ${GENE_ID}
Qry length (${QRY_LEN}) is different than Ortho length (${ORTHO_LEN}) "
fi
fa_pos.pl M441-5.chr.fa  ${ORTHO_CHR} ${ORTHO_START} ${ORTHO_END}  > ${GENE_ID}.gene.ortho_M441.fa
cat ${GENE_ID}.gene.fa ${GENE_ID}.gene.ortho_M441.fa >${GENE_ID}.combine.fa
muscle -in ${GENE_ID}.combine.fa  -clw  > ${GENE_ID}.combine.fa.aln
