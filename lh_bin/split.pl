#!/usr/bin/perl -w

if(@ARGV <2)
{
	print "$0 file split_count \n";
	exit;
}

$filename=$ARGV[0];
$split_count = $ARGV[1];

open IN,$filename or die;

while(<IN>)
{
	$line_count++;
}

close IN;
open IN,$filename or die;

my $line_index=0;
my $file_index = 0;
$every_lines = int($line_count/$split_count) + 1;
while(<IN>)
{
	if($line_index % $every_lines == 0)
	{
		open OUT, ">$filename"."_$file_index" or die;
		$file_index++;
	}
	print OUT $_;
	$line_index++;
}

