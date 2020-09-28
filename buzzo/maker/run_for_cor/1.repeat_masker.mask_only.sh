#!/bin/bash

set -euxo  pipefail

species=viridiplantae
species=cucumber
latin="Oryza sativa"
THREADS=100

mkdir -p step1_repeatmask
cat input/ref/*.fa >input/ref.fa
#cat input/peps/*.fa >input/peps.fa
#cat input/transcripts/*.fa >input/transcripts.fa

export PATH="~/lh/anaconda3/envs/repeatmodeler/bin":$PATH
export PERL5LIB=~/lh/anaconda3/envs/repeatmodeler/lib/site_perl/5.26.2/

cd step1_repeatmask
ln -s ../input/ref.fa


export PATH="~/lh/anaconda3/envs/repeatmodeler/bin":$PATH
##Remove after release

#一定要先建立目录
#One need to add Repeat Modeler to path
#
mv custom_lib.out custom_lib.out.bak
mkdir -p custom_lib.out
export PATH=~/lh/anaconda3/envs/repeatmodeler/bin/:$PATH
export PERL5LIB=""

bsub -q Q104C512G_X4  -o output.%J -e error.%J "RepeatMasker -lib ref-families.fa ref.fa -pa $THREADS -dir custom_lib.out" |awk '{print $2}' | sed 's/>//;s/<//' >>job.id

mv species_lib.out species_lib.out.bak
mkdir -p species_lib.out
bsub -q Q104C512G_X4  -o output.%J -e error.%J "RepeatMasker ref.fa -species $species -pa $THREADS -dir species_lib.out" |awk '{print $2}' | sed 's/>//;s/<//' >>job.id

