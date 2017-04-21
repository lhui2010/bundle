#35	cusau_C000N	9630428	9732199	iponi_chr	1	27497	497837
#454	cusau_C010N	2084816	2659136	iponi_chr	1	1760498	396830
#39	cusau_C000N	9105640	9419085	iponi_chr	1	3518661	679617
#742	cusau_C021N	697774	1603014	iponi_chr	1	3555566	6226657
#42	cusau_C000N	8993777	9157350	iponi_chr	1	6818176	4053928
#887	cusau_C027N	1038910	2434993	iponi_chr	1	7214421	3630309
#40	cusau_C000N	10062840	10130617	iponi_chr	1	8388319	8062936
#1003	cusau_C034N	319848	1051311	iponi_chr	1	11291082	12974411
#31	cusau_C000N	9854258	9950158	iponi_chr	1	12019851	12822717
#32	cusau_C000N	9962139	10055527	iponi_chr	1	14037601	14799757

while(<>)
{
    my @e=split;
    $hash{$e[1]}{$e[5]}  += abs($e[3] - $e[2]);
}


print "ContigID";
print "\tchr$_" for (1..15);
print "\tPutativeChr\n";
for my $k (sort keys %hash)
{
    print $k;

    my $max = 0;
    my $max_chr;
    for my $j (1..15) #sort keys %{$hash{$k}})
    {
        unless( exists ($hash{$k}{$j}))
        {
            print "\t0";
            next;
        }
        if ($max < $hash{$k}{$j})
        {
            $max = $hash{$k}{$j};
            $max_chr = $j;
        }
        print "\t",$hash{$k}{$j};
    }
    print "\t$max_chr\n";
}
