#!/usr/bin/perl -w

open GFF, "/home/SCE/liuhui/data/rice/genome/GFF3_representative/build5_with_predicted_gene.gff3" or die;

exit if (@ARGV <1);

if(@ARGV == 3)
{
	$chr = shift @ARGV;
	$start = shift @ARGV;
	$end = shift @ARGV;
	push @data, [$chr, $start, $end];
}
else
{
	while(<>)
	{
		chomp;
		@e=split;
		push @data, [@e];
	}
}


#chromosome03    build5_rep      mRNA    2806614 2808537

while(<GFF>)
{
	next unless (/mRNA/);
	@e=split;

	for my $loci(@data)
	{
		my ($chr, $start, $end) = @$loci;
		print $_ if ($e[3] >=$start and $e[4] <= $end and $e[0] eq $chr);
	}
}
