#!/usr/bin/env perl
#use strict;

while(<>)
{
	chomp;
	@e=split/\t/, $_;
#	next unless ($e[2] eq 'mRNA' or $e[2] eq 'transcript');
	next if(/^#/);
	/ID=(.*?);/;
	$name = $1;
    $name =~s/.*://;
    #1_131288720_131546951
#if gff is a cut out region
    if($e[0]=~/_/)
    {
        @list = split/_/, $e[0];
        $offset = $list[1] -1;
        $e[3] += $offset;
        $e[4] += $offset;
        $e[0] = $list[0];
    }

    print $e[0];
	for $e( $e[3], $e[4], $name, 0, $e[6])
	{
		print "\t",$e;
	}
	print "\n";
}
