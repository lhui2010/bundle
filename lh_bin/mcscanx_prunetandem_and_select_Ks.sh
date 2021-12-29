#!/bin/bash

# set -euxo pipefail

for g in *.gff;
do
    i=${g%.gff}
    detect_collinear_tandem_arrays -g ${i}.gff -b ${i}.blast -c ${i}.collinearity -o ${i}.tdarray
    perl prune_tandem.pl ${i}.tdarray ${i}.collinearity > ${i}.notd.ortho 
    selectItemH.pl ${i}.notd.ortho ${i}.ortho.kaks > ${i}.notd.ortho.kaks 
    get_ks_peak.py ${i}.notd.ortho.kaks > ${i}.notd.ortho.kaks.GMMpeak
done
