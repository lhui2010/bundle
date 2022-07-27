my %age_hash = (
    "N0"  =>   12,  
    "N1"  =>   11, 
    "N2"  =>   10, 
    "N3"  =>   9, 
    "N5"  =>   8,
    "N8"  =>   7,
    "N12" =>   6,
    "N13" =>   5,
    "N15" =>   4,
    "N18" =>   3,
    "N22" =>   2,
    "N24" =>   1);

while(<>)
{
    s/\n//;
    my $gene_age;
    #if(/\t(N\d+)\t/)
    my @e=split/\t/, $_;
    my $node_label = $e[3];
    if(exists $age_hash{$node_label})
    {
        $gene_age = $age_hash{$node_label}
    }
    else
    {
        $node_label =~s/N//;
        if($node_label > 24 or $node_label == "")
        {
            $gene_age = "0";
        }
        else
        {
            $gene_age = "NA";
        }
    }
    #$_ .= "\t\'".$gene_age."\n";
    $e[3] = $gene_age;
    $_= join("\t", @e)."\n";
    #$_ .= "\t\'".$gene_age."\n";
    print ;
}
