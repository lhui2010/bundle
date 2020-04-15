#!/usr/bin/env perl
use Bio::SeqIO;
use warnings;

if(@ARGV<1)
{
	        exit;
}

my $aa = Bio::SeqIO->new(-file=> "$ARGV[0]",
                -format=> 'Fasta' );
my %length;
while (my $seq = $aa->next_seq())
{
	        $length{$seq->id} = $seq->length;
}

print $_,"\t",$length{$_},"\n" for (sort keys %length);

