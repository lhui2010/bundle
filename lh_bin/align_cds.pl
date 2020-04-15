#!/usr/bin/env perl

use warnings;
use strict;

my @genes = @ARGV;

my $cds="/lustre/home/liuhui/project/buzzo/OrthoFinder/complete/maize_v1.1/cds/total.cds";

for my $i (@genes)
{
    system("select_fastav1.pl $i $cds >>align.fa");
}

system ("muscle -in align.fa -clw -out align.fa.clw")
