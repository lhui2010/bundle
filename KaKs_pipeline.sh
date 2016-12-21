#!/bin/bash

ORTHO=ip_so.Soly
CDS=solly.cds
PEP=solly.pep

ParaAT.pl -h $ORTHO -n $CDS -a $PEP -p proc -o ParaAT.out -f axt -k
