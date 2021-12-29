
my %self;

# Calculated as (wgda + wgdb)/2 - divab, positive indicate common
my %commonWGD;
open FH, $ARGV[0] or die;
while(<FH>)
{
    chomp;
    my @e=split;
    my ($qry, $ref) = split/\./, $e[0];
    if($qry eq $ref)
    {
        $self{$qry} = $e[1];
    }
}
close FH;

open FH, $ARGV[0] or die;
while(<FH>)
{
    chomp;
    my @e=split;
    my ($qry, $ref) = split/\./, $e[0];
    if($qry ne $ref)
    {
        $commonWGD{"$qry\t$ref"} = ($self{$qry} + $self{$ref})/2 - $e[1];
    }
}
close FH;

for my $k (sort keys %commonWGD)
{
    print("$k\t$commonWGD{$k}\n");
}
