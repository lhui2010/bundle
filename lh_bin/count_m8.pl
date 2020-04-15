#!/usr/bin/perl -w
my (%query);
while(<>)
{
	@e=split;
	$query{$e[0]}=1;
}

for $_ ( sort keys %query)
{
	$count++;
}

print $count,"\n";
