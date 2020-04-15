#!/usr/bin/perl -w
open Shorten, ">$ARGV[0].shorten" or die;
open Sdetail, ">$ARGV[0].shorten.detail" or die;
open Enlongation, ">$ARGV[0].enlongation" or die;
open Edetail, ">$ARGV[0].enlongation.detail" or die;

open IN, $ARGV[0] or die;
while(<IN>)
{
	next if (/same/);
	chomp; my @e=split;

	next if (/U<->U;/);
#ref 'gene died in resequenced rice
	if( /.<->U;/ )
	{
#		print Shorten $e[6], "\t", $e[8], "\t", $e[10], "\n";
		$isshorten{$e[6]} = 1;
		$detail{$e[6]} = $_;
	#	print $e[6], "\t", $e[8], "\t", $e[10], "\n";
	}
}

close IN;
open IN, $ARGV[0] or die;

while(<IN>)
{
	next if (/same/);next if (/U<->U;/);
	chomp; my @e=split;

	if( /U<->.;/ and not exists $isshorten{$e[6]})#once dead , no way come to life again
	{
		$isenlongation{$e[6]} = 1;
		$detail{$e[6]} = $_;
	}
}

for $key(sort keys %isshorten)
{
	print Shorten $key,"\n";
	print Sdetail $detail{$key}, "\n";
}

for $key(sort keys %isenlongation)
{
	print Enlongation $key, "\n";
	print Edetail $detail{$key}, "\n";
}
	
