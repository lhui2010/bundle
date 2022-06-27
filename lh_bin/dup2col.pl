

my $qry = $ARGV[0];
my %hash;
$_=<>;
chomp;
my @f=split;
#my $qry = $f[1];
while(<>)
{
    chomp;
    my @e=split;
    if(!exists($hash{$e[0]}))
    {
        $hash{$e[0]} = $e[1];
    }
    else
    {
        $hash{$e[0]} .= ",".$e[1];
    }
}

for my $k (sort keys %hash)
{
    print "$k\t$qry\t$hash{$k}\n";
}



