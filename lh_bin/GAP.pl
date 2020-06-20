#!/usr/bin/perl -w
use strict;
#a simple script to count gap's percentage in a multi-fasta file

my($sum, $sumN,  %length, %gap_ratio);



$/ = '>';
my $name = $ARGV[0];



while(<>)
{
        chomp;

        my ($head, $seq) = split /\n/, $_, 2;

        next unless($head && $seq);

        $seq =~ s/\n//g;

        my $r=()= $seq =~ /N/g;

#		my $r2=()=$seq =~ /[atcgn]/g;
		my $r2=()=$seq =~ /n/g;

        $sumN+=$r+$r2;
        $sum+=length($seq);

		
		#print STDERR int(100*$./$total),"%\n\b";
        #print $head, "\t", int($r*100/length($seq)), "\n";
}

print "$name\tN%:\t", int($sumN*100/$sum), "%\n$sumN\t$sum\n";

