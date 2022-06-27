#!/usr/bin/env perl
#
#Input EG. right column is the reference column
#==> Castanospermum_australe.Cercis_chinensis.ortho <==
#CASAUS_g31283.t1_Caaus    CECHI00001252-t1_Cechi
#CASAUS_g31287.t2_Caaus    CECHI00001259-t1_Cechi
#CASAUS_g31303.t1_Caaus    CECHI00001266-t1_Cechi
#CASAUS_g31307.t1_Caaus    CECHI00001274-t1_Cechi
#CASAUS_g31317.t1_Caaus    CECHI00001279-t1_Cechi
#CASAUS_g31332.t2_Caaus    CECHI00001295-t1_Cechi
#CASAUS_g31339.t1_Caaus    CECHI00001302-t1_Cechi
#CASAUS_g31342.t1_Caaus    CECHI00001307-t1_Cechi
#CASAUS_g31347.t1_Caaus    CECHI00001309-t1_Cechi
#CASAUS_g31348.t1_Caaus    CECHI00001311-t1_Cechi

my $qry = $ARGV[0];
my %hash;
$_=<>;
chomp;
my @f=split;
#my $qry = $f[1];
my %count;
while(<>)
{
    chomp;
    my @e=split;
    next if ($e[1] =~/,/);
    if(!exists($hash{$e[1]}))
    {
        $hash{$e[1]} = $e[0];
        $count{$e[1]} ++;
    }
    else
    {
        $hash{$e[1]} .= ",".$e[0];
        $count{$e[1]} ++;
    }
}

for my $k (sort keys %hash)
{
    print "$k\t$hash{$k}\n" if $count{$k} == 2;
}



