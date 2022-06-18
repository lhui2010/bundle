#!/usr/bin/env perl
use strict;
use warnings;

open DICT, $ARGV[0] or die;
open IN, $ARGV[1] or die;

my %hash;

$_=<DICT>;
#query_name
chomp;
my @e=split/\t/, $_,2;
$hash{$e[0]} = $e[1];
#$hash{"GeneID"} = $e[1];

$NA_string = "";
while(<DICT>)
{
    chomp;
    my @e=split/\t/, $_,2;
#remove transcripts
    #$e[0]=~s/_T\d+$//;
    $hash{$e[0]} = $e[1];

	if($NA_string eq "")
	{
		my @f=split/\t/, $e[1];
		$NA_string = "\tNA" * int(@f);
	}
}

while(<IN>)
{
    chomp;
    my @e=split/\t/, $_, 2;
    $e[0]=~s/^gene://;

    $_.="\t$hash{$e[0]}" if exists ($hash{$e[0]});

    $_.="$NA_string\n";

    print;
}
