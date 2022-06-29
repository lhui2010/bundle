#!/usr/bin/env perl
#
#Input EG. right column is the reference column
#Orthogroup      Cercis_chinensis        Abrus_precatorius
#OG0000001       CECHI00020547-t1_Cechi  rna-XM_027492434.1_Abpre, rna-XM_027495553.1_Abpre
#OG0000001       CECHI00018422-t1_Cechi, CECHI00018029-t1_Cechi  rna-XM_027501663.1_Abpre, rna-XM_027501689.1_Abpre
#OG0000001       CECHI00020057-t1_Cechi  rna-XM_027512481.1_Abpre, rna-XM_027512949.1_Abpre, rna-XM_027513021.1_Abpre, rna-XM_027489934.1_Abpre, rna-XM_027482171.1_Abpre, rna-XM
#OG0000003       CECHI00011898-t3_Cechi  rna-XM_027481729.1_Abpre, rna-XM_027475488.1_Abpre
#OG0000003       CECHI00017304-t1_Cechi, CECHI00017303-t1_Cechi, CECHI00010549-t1_Cechi  rna-XM_027492409.1_Abpre, rna-XM_027487539.1_Abpre
#OG0000003       CECHI00010548-t1_Cechi, CECHI00010547-t1_Cechi  rna-XM_027492408.1_Abpre, rna-XM_027484157.1_Abpre
#OG0000003       CECHI00010546-t1_Cechi  rna-XM_027492406.1_Abpre, rna-XM_027492407.1_Abpre
#OG0000003       CECHI00013641-t1_Cechi  rna-XM_027485131.1_Abpre
#OG0000003       CECHI00013261-t1_Cechi  rna-XM_027507901.1_Abpre, rna-XM_027507389.1_Abpre
#OG0000003       CECHI00013260-t1_Cechi  rna-XM_027508759.1_Abpre, rna-XM_027493744.1_Abpre
#OG0000003       CECHI00009903-t1_Cechi  rna-XM_027494499.1_Abpre, rna-XM_027478234.1_Abpre
#OG0000003       CECHI00010634-t1_Cechi  rna-XM_027487692.1_Abpre
#OG0000003       CECHI00002886-t1_Cechi  rna-XM_027512995.1_Abpre
#
my $count_per_group =2;
if (int(@ARGV)>1)
{
    $count_per_group = pop @ARGV;
}


my $qry = $ARGV[0];
my %hash;
$_=<>;
chomp;
my @f=split;
#my $qry = $f[1];
my %count;
$_=<>;
while(<>)
{
    s/ //g;
    chomp;
    my @e=split/\t/, $_;
    next if ($e[1] =~/,/);

    my @qry_list = split/,/, $e[2];
    if(int(@qry_list) == $count_per_group)
    {
        $hash{$e[1]} = $e[2];
        $count{$e[1]} +=$count_per_group;
    }
}

for my $k (sort keys %hash)
{
    print "$k\t$hash{$k}\n" if $count{$k} == $count_per_group;
}



