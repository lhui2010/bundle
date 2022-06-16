#!/usr/bin/env perl
#
$max_gene_per_group = 10;

while(<>)
{
    chomp;
    #my @e=split/\t/, $_;
    s/ //g;
    my @e = split/\t/, $_;
    if(int(@e)>1)
    {
        my @e1 = split/,/,$e[0];
        my @e2 = split/,/,$e[1];
        next if (int(@e1)<1 or int(@e2)<1 );

        for my $e1 (0..int(@e1)-1)
        {
            for my $e2(0..int(@e2)-1)
            {
                print($e1[$e1],"\t",$e2[$e2],"\n");
            }
        }
        
    }
    else
    {
        my @e1 = split/,/,$_;
        next if (int(@e1)<=1);

        for my $e1 (0..int(@e1)-1)
        {
            my $f=$e1+1;
            for my $e2($f..int(@e1)-1)
            {
                print($e1[$e1],"\t",$e1[$e2],"\n");
            }
        }
    }
}
