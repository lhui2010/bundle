#!/usr/bin/perl -w

my $count = -1;
my $line = 0;
my @tigid;
my @list;

my $file_len = shift;

my %len;
open IN, $file_len;
while(<IN>)
{
    chomp;
    my @e=split;
    $len{$e[0]} = $e[1];
}


my %strand;
$strand{'-'} = -1;
$strand{'+'} = 1;


while(<>)
{
    chomp;
    if(/::/)
    {
        $count ++;
        $_=<>;
        $_=<>;
        $_=<>;
    }
    
    $line ++;
    my @e=split;
    push @tigid, $e[0];

    my $tmp_id = $line * $strand{$e[1]};
#    print ($line, "\n", $e[1], '\n', $strand{$e[1]}, "\n", $tmp_id, "\n");

    $list[$count] .= $tmp_id. " ";

#    >tig00002059|arrow_np1212 46 12648
#    -38 -31 28 -3 -17 25 -7 -8 12 15 -11 -2 35 -41 -36 -30 39 21 23 -34 26 29 33 -27 37 -42 -22 -24 32 20 18 19 1 5 10 16 14 9 -13 6 4
#    40

}

for my $i (0..$#tigid)
{
    print ">", $tigid[$i], "\t", $i + 1, "\t", $len{$tigid[$i]}, "\n";
}

for my $c (@list)
{
    print $c, "\n";
}
