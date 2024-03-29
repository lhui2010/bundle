#!/usr/bin/env perl


#cu00000181      .       mRNA    264475  272937  .       +       .       ID=evm.model.tig00000181.16.2;Parent=evm.TU.tig00000181.16;Name=EVM%20prediction%20tig00000181.16

while(<>)
{
    chomp;
    my @e=split;
    next unless ($e[2] eq mRNA);
    chomp;
    if(/ID=(.*?);/ or /ID=(.*)$/ )
    {
        $mRNAID = $1;
        $start{$mRNAID} = $e[3];
        $end{$mRNAID} = $e[4];
        $ctg{$mRNAID} = $e[0];
    }
    else
    {
        next;
    }
}

for my $k (sort keys %start)
{
    print $ctg{$k},"\t",$k,"\t", $start{$k}, "\t", $end{$k},"\n";
}
    
