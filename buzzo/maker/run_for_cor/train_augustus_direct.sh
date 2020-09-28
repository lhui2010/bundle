#!/bin/bash
#BSUB -J direct_train_aug      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail 

ROOT=$PWD


WORKING_DIR=${ROOT}
NUMSPLIT=500

#coriariav2.0 is augustus directly predicted, while coriaria_v2.0 is busco predicted
AUGUSTUS_SPECIES_NAME=coriariav2.0_direct

CDNA_FASTA=cdna.fasta

export PATH=~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/:$PATH



NUMFOUND="`grep -c '>' uni.ann`"
if [ ${NUMFOUND} -gt 499 ]; then
    NUMFOUND=500
fi


TEMPSPLIT=$((NUMFOUND/2))
NUMSPLIT=${TEMPSPLIT/.*}
echo ${NUMSPLIT}
perl GC_specific_MAKER/fathom_to_genbank.pl --annotation_file uni.ann  --dna_file uni.dna  --genbank_file augustus.gb  --number ${NUMFOUND}
perl -e  'while (my $line = <>){ if ($line =~ /^LOCUS\s+(\S+)/) { print "$1\n"; } }'  ${WORKING_DIR}/augustus.gb  >  ${WORKING_DIR}/genbank_gene_list.txt 
perl GC_specific_MAKER/get_subset_of_fastas.pl   -l  ${WORKING_DIR}/genbank_gene_list.txt    -f ${WORKING_DIR}/uni.dna  -o  ${WORKING_DIR}/genbank_gene_seqs.fasta

perl ~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/randomSplit.pl ${WORKING_DIR}/augustus.gb ${NUMSPLIT}


perl ~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/autoAug.pl --species=$AUGUSTUS_SPECIES_NAME --genome=${WORKING_DIR}/genbank_gene_seqs.fasta --trainingset=${WORKING_DIR}/augustus.gb.train --cdna=$CDNA_FASTA  --noutr

cd ${WORKING_DIR}/autoAug/autoAugPred_abinitio/shells

# The number of ./aug# scripts is variable.  Run until the new file does not exist.
x=1
while [ -e ./aug${x} ]
do
    echo "A.  $x"
    ./aug${x} &
    let x=x+1
done

wait

cd $WORKING_DIR

perl ~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/autoAug.pl --species=$AUGUSTUS_SPECIES_NAME --genome=${WORKING_DIR}/genbank_gene_seqs.fasta --useexisting --hints=${WORKING_DIR}/autoAug/hints/hints.E.gff  -v -v -v  --index=1

cd ${WORKING_DIR}/autoAug/autoAugPred_hints/shells/

let x=1
while [ -e ./aug${x} ]
do
    echo "B.  $x"
    ./aug${x} &
    let x=x+1
done

wait

cd ${WORKING_DIR}

augustus --species=$AUGUSTUS_SPECIES_NAME augustus.gb.test
