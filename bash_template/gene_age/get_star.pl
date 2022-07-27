
#my %star_grouped;
my %star_grouped;
my %star_tax_covered;
my %star_contain_dup;

my $total_taxa_number = 47;
my $taxa_number_cutoff = 23;
my $dup_cutoff = 2;

while(<>)
{
    chomp;

    my @e=split;
    $star_grouped{$e[0]} = 1;
    if($e[-1] ne "NA")
    {
        $star_tax_covered{$e[0]} ++;
    }

    if($e[-1] eq "duplicate")
    {
        $star_contain_dup{$e[0]} ++;
    }
}

for my $k(sort keys %star_grouped)
{
    print $k, "\t";
    if (exists $star_grouped{$k})
    {
        print "@-@";
    }
    print "\t";
    if ($star_tax_covered{$k} > $taxa_number_cutoff)
    {
        print "@-@";
    }
    print "\t";
    if($star_contain_dup{$k} > $dup_cutoff)
    {
        print "@-@";
    }
    print "\n";
}

