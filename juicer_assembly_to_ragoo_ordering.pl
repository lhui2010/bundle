
#>000263F|arrow_np1212 268 2550
#>000264F|arrow_np1212 269 2356
#58 -90 -27 198 -95 75 -144 35 42
#
#perl -w juicer_assembly_to_ragoo_ordering.pl falcon_v340_sgs_polish.1.review.assembly
my $chr_id = 1;
my %chr_table;
while(<>)
{
    chomp;
    if(/^>/)
    {
        my @e=split;
        $e[0]=~s/>//;
        $index2contig_name{$e[1]} = $e[0];
    }
    else
    {
        my $chr_key = 'chr'.sprintf("%02d", $chr_id);
        #print $chr_key; exit;
        my @e=split;

#tig00000809|arrow_np1212        +       1.0     1.0
#tig00000814|arrow_np1212        -       1.0     1.0
#tig00000812|arrow_np1212        +       1.0     1.0
#tig00000807|arrow_np1212        -       0.9949528219634588      1.0
        for my $e(@e)
        {
            my $strand = "+";
            if($e < 0)
            {
                $strand = "-";
                $e = 0 - $e;
            }
            $chr_table{$chr_key}.=join("\t", ($index2contig_name{$e}, $strand, '1.0', '1.0'))."\n";
        }
        $chr_id ++;
    }
}

`mkdir -p ragoo_ordering`;
for my $i (sort keys %chr_table)
{
    open OUT, ">ragoo_ordering/".$i."_orderings.txt" or die;
    print OUT $chr_table{$i};
    close OUT;
}

