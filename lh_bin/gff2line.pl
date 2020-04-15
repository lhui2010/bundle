#!/usr/bin/perl -w
while(<>)
{
	next unless (/mRNA/);
	/ID=(.*)-..;Name.*Note=(.*);Transcript/;
	if (defined $2)
	{
		$gff{$1} = $2;
	}
	else
	{
		/ID=(.*);Name.*Note=(.*) n=/;
		print $1, "\t", $2, "\n";
	}
}

for $key (sort keys %gff)
{
	print $key, "\t", $gff{$key}, "\n";
}
