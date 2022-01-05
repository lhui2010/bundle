#!/usr/bin/perl -w


die "sumdict.pl table(two column with first as key and second as value) > table's sum" if @ARGV < 2;

#Input
#abc 1
#abc 2

#Output
#abc 3

my %sum;

while(<>)
{
	my @e=split;
	$sum{$e[0]}+=$e[1];
}

for my $k (sort keys %sum)
{
    print $k, "\t", $sum{$k}, "\n";
}
