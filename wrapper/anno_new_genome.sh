#!/bin/bash

set -eo

# User defined variables
PREFIX=CASAUS
REF=Castanospermum_australe.fa
PEP=arath_med_sprot.pep
ISOSEQDIR=workdir_isoseq_clpla
# User defined variables end

python -m iga.annotation.repeat repeatmasker --species Viridiplantae --denovo F ${REF} &
python -m iga.annotation.repeat repeatmasker --species '' --denovo T ${REF} &


wait
python -m iga.annotation.repeat post_repeatmasker workdir_repeatmask_${REF} ${REF}

ln -s $PWD/workdir_repeatmask_${REF}/Full_mask/full_mask.complex.reformat.gff3 ${REF}.repeat.gff

REPEAT_GFF=${REF}.repeat.gff

MASKEDREF=${REF}.masked.fa

if [ ! -e  ${MASKEDREF} ]
then
    echo "Error, masked genome not exists, check log"
    exit
fi

### RUN braker

source activate braker

export GENEMARK_PATH=/nfs/liuhui/bin/braker/dependency/gmes_linux_64
export PROTHINT_PATH=/nfs/liuhui/bin/braker/dependency/ProtHint/bin
export CDBTOOLS_PATH=/nfs/liuhui/bin/braker/dependency/cdbfasta
export PATH=/nfs/liuhui/bin/braker/dependency/GUSHR:$PATH
export MAKEHUB_PATH=/nfs/liuhui/bin/braker/dependency/MakeHub-1.0.5
export PATH=/nfs/liuhui/bin/braker/dependency:$PATH
export PATH=/nfs/liuhui/bin/braker/dependency/Augustus/bin:/nfs/liuhui/bin/braker/dependency/Augustus/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=/nfs/liuhui/bin/braker/dependency/Augustus/config
export PATH=/nfs/liuhui/bin/braker/BRAKER/scripts:$PATH

module load boost

PEP=odb10_plants.fasta
THREADS=48

for GENOME in $MASKEDREF
do
    while [ ! -e ${GENOME} ]
    do
        echo "${GENOME} is not yet generated"
        sleep 1m
    done
    bsub  -R "span[hosts=1]" -q Q104C512G_X4  -o output.%J -e error.%J -J braker -n 48 "braker.pl -gff3 --cores=${THREADS} --genome=${GENOME} --prot_seq=$PEP --softmasking --workingdir braker_${GENOME} && pushd braker_${GENOME} && longest_transcript_from_gff.py braker.gff3 && sed "s/file_1_file_1/${PREFIX}/g" braker.gff3.longest > ${REF%.fa}.gff3 && cp ../${MASKEDREF} ${REF} && add_suffix ${REF} ${REF%.fa}.gff3 && for i in *.new; do mv ${i} ${i%.new}; done && popd "
done


