#!/usr/bin/perl

my $usage = "usage: $0 gene.bed repeat.bed gene_intersect_repeat.bed\n\n";

my $overlap_bed = pop @ARGV or die $usage;
open INTEC_BED, $overlap_bed or die;

my %fa_len;
#xx.bed bb.bed
while(<>)
{
	my @e=split;
    unshift @e, "";
    my $key = $e[4];
    my $len = $e[3] - $e[2];
    $fa_len{$key} += $len;
}
    

my %sum_overlap;
while(<INTEC_BED>)
{
#	/size([0-9]*)/;
# 8 9 10 13
	my @e=split;
    unshift @e, "";
    $key = $e[10]."\t".$e[4];
	$sum_overlap{$e[10]}{$e[4]}+=$e[13];
}

for my $k (sort keys %sum_overlap)
{
    for my $t (sort keys %{$sum_overlap{$k}})
    {
        print STDERR join("\n", ($k, $t, $fa_len{$k}, $fa_len{$t})) if $fa_len{$k} * $fa_len{$t} ==0;
        print join("\t", ($k, $t, $sum_overlap{$k}{$t}, $fa_len{$k}, $fa_len{$t}, $sum_overlap{$k}{$t} / $fa_len{$k}, $sum_overlap{$k}{$t} / $fa_len{$t})), "\n";
    }
}
#    print $k, "\t", $sum_overlap{$k}, "\t", $sum_cds{$k}, "\t", $sum_orig_cds{$k}, "\t", $sum_overlap{$k}/($sum_cds{$k}+1), "\t", $sum_overlap{$k}/($sum_orig_cds{$k}+1), "\n";
#    print STDERR $k if $sum_cds{$k} * $sum_overlap{$k} ==0;
   # print $k, "\t", $sum_overlap{$k}, "\t", $sum_cds{$k}, "\t", $sum_overlap{$k}/$sum_cds{$k}, "\n";
