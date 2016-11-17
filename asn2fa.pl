

#note="scaffold: scaffold0001"
#ORIGIN
#        1 ccctcgtgcc aaggccgccg cgacgattcg cgacgaccgt ggaacgtggc gatgccttct

our ($seq, $mark, %fa);

while(<>)
{
    if(/scaffold: (.*)"/)
    {
        $scaff_id=$1;
    }
    if(/ORIGIN/)
    {
        &init_read;
        next;
    }
    if(/^LOCUS/)
    {
        &end_read;
    }
    if($mark)
    {
        &read_seq;
    }
}

for my $k (sort keys %fa)
{
    next if ($k eq "");
    print ">",$k,"\n";
    print(uc($fa{$k}),"\n");
}


sub init_read
{
    $seq= "";
    $mark = 1;
}

sub end_read
{
    $fa{$scaff_id} = $seq;
    $mark = 0;
}

sub read_seq
{
    chomp;
    s/\s+//g;
    s/^\d+//;
    $seq.=$_;
}

