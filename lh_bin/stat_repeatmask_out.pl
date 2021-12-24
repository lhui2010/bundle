

#  4942    7.3  2.7  1.1  iponi_chr1                 1      733 (42878941) C rnd-1_family-381   Unknown             (1696)    745       1       1  

$_=<>;
$_=<>;
$_=<>;

while(<>)
{
    chomp;
    my @e=split;
    
    $hash{$e[10]}+= $e[6]-$e[5];
}

$genome_size=734803496;

for my $k (sort {$hash{$b}<=>$hash{$a}} keys %hash)
{
    print $k, "\t\t\t",  $hash{$k}/$genome_size, "\t", $hash{$k}, "\n";
}

