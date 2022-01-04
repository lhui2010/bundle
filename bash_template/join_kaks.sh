#!/bin/bash

cd mcscanX
for i in `find . -mindepth 1 -maxdepth 1 -type d -name '*.ParaAT.out'`
do
    echo -e "Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites(0:2:4)\tSubstitutions\tS-SubstitutionsN-Substitutions\tFold-S-Substitutions(0:2:4)\tFold-N-Substitutions(0:2:4)\tDivergence-Time\tSubstitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\tGC(1:2:3)\tML-Score\tAICc\tAkaike-Weight\tModel" > ${i%.ParaAT.out}.kaks
    find ${i} -name "*.kaks" |xargs tail -n +2  >> ${i%.ParaAT.out}.kaks
done

