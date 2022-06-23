#!/usr/bin/env perl


#0.7
my $single_cutoff=34;
#0.5 cutoff=25
if(int(@ARGV) > 1)
{
    $single_cutoff = pop;
}

#header
$_=<>;
my @e=split/\t/,$_;
my $tag = shift @e;
shift@e;
shift@e;
#unshift @e, "SingleCount";
unshift @e, "HOG";
print(join("\t", @e));

while(<>)
{
    chomp;
    my @result_genes;
    my $single_count = 0;
    my @e=split/\t/; 
    my $hog = shift @e; 
    shift @e; 
    shift @e; 

    for my $e(@e)
    {
        unless($e=~/,/ )
        {
            if($e ne "")
            {
                $single_count ++;
            }
            else
            {
                $e=" ";
            }
            push @result_genes, $e;
        }
        else
        {
            push @result_genes, "";
        }

    } 
    if ($single_count >= $single_cutoff)
    {
        #print "$hog\t$single_count\t", join("\t", @result_genes), "\n";
        print "$hog\t", join("\t", @result_genes), "\n";
    }
    else
    {
        #print STDERR $single_count;
    }
}
