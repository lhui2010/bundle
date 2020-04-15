#!/usr/bin/perl -w

open IN, $ARGV[0] or die;
shift @ARGV;
for my $i(@ARGV)
{
	$line_number{$i} = 1;
}

while(<IN>)
{
	if (exists $line_number{$.})
	{
		print;
	}
}
