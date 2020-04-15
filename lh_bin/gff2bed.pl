#!/usr/bin/env perl
#use strict;

while(<>)
{
	chomp;
	@e=split/\t/, $_;
	next unless ($e[2] eq 'mRNA' or $e[2] eq 'transcript');
	/ID=(.*?);/;
	$name = $1;
    $name =~s/.*://;
    print $e[0];
    $e[3]-=1;
	for $e( $e[3], $e[4], $name, 0, $e[6])
	{
		print "\t",$e;
	}
	print "\n";
}
