#!/usr/bin/env perl

while(<>)
{
    next unless (/rank:1/);
my @e=split/\|/, $_;
my ($chr, $loci) = split/:/, $e[1];
my ($start, $end) = split/\.\./, $loci;
my $gene_name = $e[0];
my $strand = $e[2];
my $score = $e[4];
$score=~s/score://;

print "$chr\t$start\t$end\t$gene_name\t".int($score)."\t$strand\n";
}
