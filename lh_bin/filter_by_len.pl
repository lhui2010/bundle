#!/usr/bin/perl -w
use strict;

my($seq, $head, $r, $sum, $sumGC, $seq2, $LENGTH);

print "eg: $0 40 seq.fa >result.fa\n"if(@ARGV <2);

$LENGTH=$ARGV[0];

open IN, $ARGV[1] or die;
$/ = '>';
while(<IN>)
{
	chomp;
	
	my ($head, $seq) = split /\n/, $_, 2;

	next unless($head && $seq);

	$seq2 = $seq;
	$seq2 =~ s/\n//g;

	my $len = length ($seq2);
	
	print ">", $head, "\n", $seq if ($len >=$LENGTH);
}


