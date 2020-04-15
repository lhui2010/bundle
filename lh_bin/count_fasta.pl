#!/usr/bin/perl -w
while(<>)
{
	chomp;
	next if (/>/);
	$count += length;
}
print $count,"\n";
