#!/usr/bin/perl -w

#my $column = shift;
my $end = 3;
my $start = 2;

while(<>)
{
    chomp;
    next if ($_ eq "");
	my @e=split;
	$sum+=$e[$end] - $e[$start];
}

print "Clean sum:  ", $sum, "\n";
