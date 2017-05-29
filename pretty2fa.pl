#!/bin/perl
use strict;

#input: genewise pretty format result
#outut: pep and cds sequences in OUT/ dir

#pretty format output of genewise -> pep, cds while trailing all indels
#
#genewise output
#Score 671.53 bits over entire alignment
#Scores as bits over a synchronous coding model
#
#Warning: The bits scores is not probablistically correct for single seqs
#See WWW help for more info
#
#C000N0093G1.1      1 QSILEMANRVSSFSRLVQXVTASCLLHP-VVDAADHHDFSVVNTAGEES
#                     QSILEMANRVSSFSRLVQXVTASCLLHP VVDAADHHDFSVVNTAGEES
#                     QSILEMANRVSSFSRLVQXVTASCLLHPXVVDAADHHDFSVVNTAGEES
#C000N0093G1.1      1 caatgagaagtattccgc5gagttcccccggggggccgttggaaggggt
#                     agttatcagtcgtcgtta tcccgttactttaccaaaatcttaccgaac
#                     accgaggcgcctcctgcg tagttgtctcccccctccccccctcccaac
#

for (my $i++ < 8)
{
	my $tmp=<>;
}

my $prev = "";
my $pepID;
my $cdsID;
my $pep_seq;
my $cds_seq;

my $mark = 0;
while(<>)
{
	$mark = 1 if (/^C/);
	next unless ($mark);
	my $seq;
	chomp;
	s/\s+$//;
	if (/^C/ and $prev eq "") #start of a protein sequence
	{
		($pepID, undef, $seq) = split  /\s+/, $_, 3;
		$seq=<>;
		$seq=<>;
		$seq =~s/X//g;
		$seq =~s/\s+//g;
		$seq =~s/!//g;
		$pep_seq .=$seq;
	}
	if (/^C/ and $prev ne "") #start of a cds sequence
	{
		($cdsID, undef, $seq) = split /\s+/, $_, 3;
		$seq =~s/\d//g;
		$seq =~s/!//g;
		$seq =~s/\s+//g;
		my @codon1 = split//, $seq;
		
		$seq=<>;
		$seq =~s/\d//g;
		$seq =~s/!//g;
		$seq =~s/\s+//g;
		my @codon2 = split//, $seq;

		$seq=<>;
		$seq =~s/!//g;
		$seq =~s/\s+//g;
		my @codon3 = split//, $seq;

		for my $i (0..$#codon1)
		{
			my $codon = uc($codon1[$i].$codon2[$i].$codon3[$i]);
			$cds_seq .= $codon if($codon ne "TGA" and $codon ne "TAA" and $codon ne "TAG")# and $codon ne "ATG");
		}
	}
	$prev = $_;
}

#print protein sequences
open  OUT, ">OUT/$pepID.pep";
print OUT ">$pepID\n$pep_seq\n";
close OUT;

open  OUT, ">OUT/$cdsID.cds";
print OUT ">$cdsID\n",$cds_seq,"\n";
close OUT;
