

my $missing_cutoff=5;
my $double_cutoff=5;
my $single_cutoff=40;

#header
$_=<>;
print;

while(<>)
{
    my $missing_count = 0;
    my $double_count = 0;
    my $single_count = 0;
    my @e=split; 
    shift @e; 
    pop @e; 

    my $flag=1;
    for my $e(@e)
    {
        $flag=0 if $e > 2;
        $double_count ++ if $e == 2;
        $single_count ++ if $e == 1;
        $missing_count ++ if $e == 0;
    } 
    if ($double_count <= $double_cutoff && 
        $missing_count <= $missing_cutoff && 
        $single_count >= $single_cutoff &&
        $flag)
    {
        print;
    }
}
