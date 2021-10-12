#!/usr/bin/env perl 

# Input: *.fai

my $chr_max_number = 30;
my $chr_portion_thres = 0.85;

for my $f (@ARGV)
{
    my $sorted = `sort -k2,2nr $f`;

    my @sorted = split/\n/, $sorted;


    # The theroy is, if the 30 largest contigs/scaffolds exceed 85% of total chromosome, then this fasta is recognized as chromosomal-level assembly

    my $pseudochr_size = 0;
    my $total_size = 0;
    my $count = 0;
    
    for my $s (@sorted)
    {
        my @field = split/\s+/, $s;
        if($count < $chr_max_number)
        {
            $pseudochr_size += $field[1];
        }
        $count ++;
        $total_size += $field[1];
    }
    my $ratio = $pseudochr_size / $total_size;
    print STDERR  $f, "\t", $pseudochr_size, "\t", $total_size, "\t", $ratio, "\n";

    if ($ratio > $chr_portion_thres)
    {
        print "$f\tOK\n";
    }

}


