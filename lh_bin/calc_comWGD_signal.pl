#!/usr/bin/env perl
#
$usage = "Usage: \n$0 peak.txt \nPeakEg: Abrus_precatorius.Aeschynomene_evenia.ortho     0.6173615640569423(LeftParalogKsPeak)      0.741019588316342(RightParalogKsPeak)       0.6365870986407788(orthologKsPeak)\n Postive means common\n";
die $usage if int(@ARGV) < 1;
#print("qry\tref\tresult\n");
while(<>)
{
    next if (/^#/);
    chomp;
    my @e=split;
    my @f = split/\./, $e[0];
    print $f[0], "\t", $f[1], "\t";
    my $stat= ($e[1]+$e[2])/2 - $e[3] ;
    print $stat, "\n";
}

