#!/usr/bin/perl -w 
use strict;
#alnid    qry_contig_id         qrystart qryend  ref_chr refstart        refend   strand
#10632   tig00000000_pilon       218030  2370094 chr5    20804540        21809463 plus
#7797    tig00000000_pilon       429448  2203864 chr2    51613282        52377921 plus
#11360   tig00000000_pilon       511129  3650700 chr6    14169939        15995147 plus
#9238    tig00000000_pilon       511129  3543555 chr3    12614209        18453141 plus
#9243    tig00000000_pilon       592563  2190660 chr3    21217762        22300188 plus
#12228   tig00000000_pilon       616011  2356784 chr7    2688044         3591166  plus
#12905   tig00000000_pilon       616011  2190660 chr8    23491691        24664888 plus
#2884    tig00000000_pilon       616011  3351837 chr0    129366264       13770040 plus
#4116    tig00000000_pilon       616011  4363183 chr1    14015971        23884483 plus

$_=<>;
chomp;
my ($lid,$lqry_tig, $lqry_start, $lqry_end, $ltig, $lstart, $lend, $lstrand)  = split; #sorted
my $lprint = $_;

my %classify_chr;
my %chr;
while(<>)
{
    chomp;
    my ($id,$qry_tig, $qry_start, $qry_end, $tig, $start, $end, $strand) = split;
    
#to print header
    $chr{$tig} = 1;

    if ($qry_tig ne $lqry_tig or $tig ne $ltig or $start > $lend)
    {
        $classify_chr{$lqry_tig}{$ltig}+=$lend-$lstart;
        print $lprint,"\n";
        ($lid,$lqry_tig, $lqry_start, $lqry_end, $ltig, $lstart, $lend, $lstrand)=
        ($id,$qry_tig, $qry_start, $qry_end, $tig, $start, $end, $strand);
        $lprint = $_;
    }

    if($start < $end and $end >$lend)
    {
        $lend=$end;
    }
}

print $_,"\n";

open OUT, ">contig_chr_stat_table" or die;
print OUT "Contig";
print OUT "\t$_" for (sort keys %chr);
#also output maximum aligned chr, or putative location chr of contig
print OUT "\tmax";
print OUT "\n";
for my $k (sort keys %classify_chr)
{
    print OUT $k;
    #for my $chr (sort {$classify_chr{$k}{$b} <=> $classify_chr{$k}{$a}} keys %{$classify_chr{$k}})
    my $maxid = "chr0";
    my $maxaln= 0;
    for my $chr (sort keys %chr)
    {
        if (exists $classify_chr{$k}{$chr})
        {
            print OUT "\t$classify_chr{$k}{$chr}";
            
            if($classify_chr{$k}{$chr} > $maxaln and $chr ne "chr0") #a hack to sovle exclussive chr0 problem
            {
                $maxaln = $classify_chr{$k}{$chr};
                $maxid = $chr;
            }
        }
        else
        {
            print OUT "\t0";
        }
    }
    print OUT "\t$maxid";
    print OUT "\n";
}
