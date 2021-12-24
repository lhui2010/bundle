

my $max = 6;
my $min = 0;
my $sep = 0.1;

my %hash;

while(<>)
{
    chomp;
    my ($Sequence,$Method,$Ka,$Ks,$Ka_Ks,$P_Value,$Length,$S_Sites,$N_Sites, undef) = split /\t/, $_;
#Fold-Sites(0:2:4)	Substitutions	S-Substitutions	N-Substitutions	Fold-S-Substitutions(0:2:4)	Fold-N-Substitutions(0:2:4)	Divergence-Time	Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)	GC(1:2:3)	ML-Score	AICc	Akaike-Weight	Model
    &hist_add($Ks);
    $total_count++;
}

&hist_print();

sub hist_add()
{
    my $input = shift;
    return 0 unless ($input =~ /^\d+/);
#    $input = $max if ($input >$max);
 #   $input = $min if ($input <= $min);
    next if($input >$max or $input <=$min);
    ($input) = (int($input*10)/10);
    $hash{$input}++;
}

    

sub hist_print()
{
    for my $k (sort {$a<=>$b} keys %hash)
    {
        print $k, "\t", $hash{$k}/$total_count, "\n";
    }
}
