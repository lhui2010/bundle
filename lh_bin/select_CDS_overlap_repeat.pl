#!/usr/bin/env perl
#NW_020874290.1_Abpre    167550  169071  rna-XM_027475485.1_Abpre        0       +       1       1521    1521    1.0000000
#
#Wed Jun 29 16:39:52 CST 2022
#for i in *.gene.bed; do bedtools coverage -a $i -b ${i%.gene.bed}.repeat.gff > ${i}.xrepeat.coverage; done
#Wed Jun 29 16:40:04 CST 2022
#for i in *gff3; do bsub "python -m jcvi.formats.gff bed --type=CDS --key=Parent ${i} -o ${i%.gff3}.gene.bed"; done
#Wed Jun 29 16:41:34 CST 2022
#$for i in *.gene.bed; do bsub "bedtools coverage -a $i -b ${i%.gene.bed}.fa.repeat.gff -s > ${i}.xrepeat.coverage"; done
#
my %len_hash;
my %cov_hash;
my $cutoff = 0.5;
while(<>)
{
    my @e=split;
    my $feat_name = $e[3];
    my $feat_len = $e[-2];

    my $match_len = $e[-3];
    $cov_hash{$feat_name} += $match_len;
    $len_hash{$feat_name} += $feat_len;
}

for my $name (sort keys %cov_hash)
{
    if($cov_hash{$name} /  $len_hash{$name} > $cutoff)
    {
        print $name, "\n";
    }
}
