
#!/usr/bin/perl -w

#Aligned Fasta File to Phylip Format

$/='>';

while(<>)
{
        chomp;
        next if ($_ eq "");
        my ($info, $seq) = split /\n/, $_, 2;
        $info =~s/\s+.*//g;
        $seq =~s/\n//g;

        $hash{$info} = $seq;
        $len = length($seq);
        $count ++;
}


print "$count  $len\n";
for my $k (sort keys %hash)
{
        print "$k  $hash{$k}\n";
}

