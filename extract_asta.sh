#!/bin/bash
input=$1
grep '"1^,2^"' $input >${input}.asta_5alt.gtf
grep '"1-,2-"' $input >${input}.asta_3alt.gtf
grep '"0,1-2^"' $input >${input}.asta_exonskip.gtf
grep '"0,1^2-"' $input >${input}.asta_intron_retention.gtf
grep '"1-2^,3-4^"' $input >${input}.asta_mutually_exclusive_exons.gtf


