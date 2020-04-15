#!/usr/bin/perl -w
$p = 1e-5;

while(<>)
{
	next if (/NA/i or /NULL/i);
	chomp;
	my @e=split /\t/, $_;
	if($e[8] < $p)# or (exists $name{$e[0]} and $name{$e[0]}[4] > $e[8]))
	{
		$name{$e[0]} = join "\t", ($e[0], $e[2], $e[6], $e[7], $e[8], $e[11], $e[12], $e[13]);
	}
}
for $key (sort keys %name)
{
#	for $e( @{$name{$key}})
#	{
#		print $name{$key}[$e], "\t";
#	}
	print $name{$key};
	print "\n";
}
