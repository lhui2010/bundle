#!/usr/bin/perl -w


die "frac.pl partial_table(numerator)  total_table(denominator) > partial_table's percentage" if @ARGV < 2;


#Input
#abc 1
#abc 2

#Output
#abc 3

my %sum;

open PART, $ARGV[0] or die;
open TOTAL, $ARGV[1] or die;

my (%part, %total);

while(<TOTAL>)
{
    chomp;
	my @e=split;
	$total{$e[0]}=$e[1];
}

while(<PART>)
{
    chomp;
	my @e=split;
	$part{$e[0]}=$e[1];
    my $frac = $e[1] / $total{$e[0]};
    print $_, "\t", $frac, "\n";
}
