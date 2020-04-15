


open LEN, $ARGV[0];
open IN, $ARGV[1] or die;
while(<LEN>)
{
	chomp;
	my @e=split;
	$len{$e[0]}= $e[1];
}

while(<IN>)
{
	chomp;
	my @e=split;

	$e[0] .="\t$len{$e[0]}";

	my $print = join"\t", @e;
	print $print, "\n";
}
