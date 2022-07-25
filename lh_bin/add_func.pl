#!/usr/bin/env perl
use strict;
use warnings;

my $final_table = pop @ARGV;


#Support multiple dict
open IN, $final_table or die;

my %hash;

#open DICT, $ARGV[0] or die;
#  $_=<DICT>;
#  #query_name
#  chomp;
#  my @e=split/\t/, $_,2;
#  $hash{$e[0]} = $e[1];
#  #$hash{"GeneID"} = $e[1];

my $NA_string = "";
while(<>)
{
    chomp;
    my @e=split/\t/, $_,2;
#remove transcripts
    #$e[0]=~s/_T\d+$//;
    $hash{$e[0]} = $e[1];

	if($NA_string eq "")
	{
		my @f=split/\t/, $e[1];
		$NA_string = "\t\'" x 2;#int(@f);
	}
}

while(<IN>)
{
    chomp;
    my @e=split/\t/, $_, 2;
    $e[0]=~s/^gene://;

	if(exists ($hash{$e[0]}))
	{
    	$_.="\t$hash{$e[0]}\n";
	}
	else
	{
		$_.="$NA_string\n";
	}

    print;
}
close IN;
