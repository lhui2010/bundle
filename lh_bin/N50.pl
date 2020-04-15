#!/usr/bin/perl -w
use strict;

my ($len, $seq, %hash, $total, $adder, $half, $key, $fname);

$fname = $ARGV[0];
$/ = '>';
my $max = 0;
my $min = 10000;
while(<>)
{
	($len, $seq) = split /\n/, $_, 2;

	next unless (defined $seq and $seq ne "");

	$seq =~ s/\n//g;
	$seq =~ s/\s//g;

	$len = length($seq);

	$hash {$len} += $len;
	$total += $len;

	$max = $len if ($max <$len);
	$min = $len if ($min >$len);
}

print "spname\tmin\tmax\tN50\tN90\t";
print $fname,"\t",$min,"\t",$max,"\t";
$half = $total / 2;

my $nin_per = $total /10 * 9;

my $signal_5 = 1;
my $signal_9 = 1;


for $key(sort {$b<=>$a} keys %hash)
{
	$adder += $hash{$key};
	if( $adder >= $half and $signal_5)
	{
		print $key, "\t" ;
		$signal_5 = 0;
#		exit;
	}
	if( $adder >= $nin_per and $signal_9)
	{
		print $key, "\n";
		$signal_9 = 0;
	}
}
