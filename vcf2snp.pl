#!/bin/perl -w
use strict;

if (@ARGV < 1 )
{
    print "$0\nUsage: $0 vcf >vcf.genotype\n";
    exit;
}

my @accession = qw/W01 W03 W04 W06 W07 W13 W14 W28 W37 W38 W39 W40 W41 W42 W43 W44 W45 W46 W47 W48 W49 W50 W51 W52 W53 W54 W55 W56 W57 W58 W59 W60 W61 W62 W63 W64 W68 W69 W70 W71/;


my %het_gen = ("AC" => "M","CA" => "M","AT" => "W", "TA" => "W", "AG" => "R", "GA" => "R", "CT" => "Y", "TC" => "Y", "CG" => "S" , "GC" => "S", "GT" => "K", "TG" => "K", "AA" => "A", "TT" => "T", "CC" => "C", "GG" => "G", "-" => "-", "--" => "-");

#print header
#my $count;

while(<>)
{
	next if (/^#/);
	chomp;
	my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @tag_sp) = split;

#heterozygous site in GATK denoted as: C,A

	print $CHROM, "\t", $POS;

	for (my $i =0; $i <int(@tag_sp) ; $i++)
	{
		my $genotype= "-";
		my @alt = split (/,/, $ALT);
		unshift @alt, $REF;

#		print $alt[0],"\n", $alt[1],"\n",$alt[2],"\n"; exit;

		$genotype = $het_gen{&tag2genotype($tag_sp[$i], @alt)};

#		$genotype = &tag2genotype($tag_sp[$i],  @alt);
#		$genotype = $het_gen{&tag2genotype($tag_sp[$i], $REF, $ALT)};
		print "\t", $genotype;
	}
	print "\n";

	#exit if ($count ++>100);
}

sub tag2genotype
{
	my ($info, @geno) =  @_;
#1/1:0,14:14:42:585,42,0 ./.:.:.:.:.
	return "-" if ($info=~/^\./);

#if gatk met heterozygous site, it will mark as 2/2 instead of 1/1 for alt allel
#hah, that's only your guess and it was wrong!

	if($info =~/(\d)\/(\d)/)
	{
		my ($a, $b) = ($1, $2);
		#print $geno[$a].$geno[$b]; exit;
		my $return = $geno[$a].$geno[$b];
		return $return if (defined $return);
	}
}

