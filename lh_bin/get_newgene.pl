use strict;

open ORDER, "cu_ip.sorted.gff" or die;

#scaffold0003    INIL00g00001.t1 1       4887
#scaffold0005    INIL00g00002.t1 2       201
#scaffold0005    INIL00g00003.t1 3156    5483

$_=<ORDER>;
my @e=split;
my $last_scaf = $e[0];
my $last_gene = $e[1];

my (%next, %prev);
while(<ORDER>)
{
    chomp;
    my @e=split;

    if($e[0] eq $last_scaf)#same scaffold
    {
        $next{$last_gene} = $e[1];
        $prev{$e[1]} = $last_gene;
    }
    $last_scaf = $e[0];
    $last_gene = $e[1];
}

open SYN, "cu_ip.cross" or die;
my ($ALN_ID, %strand, %ortholog, %start, %end, %info);
while(<SYN>)
{
### Alignment 20: score=439.0 e_value=1.6e-17 N=9 scaffold0019&tig00000246_pilon plus
#INIL01g00028.t1 evm.model.tig00000246_pilon.942

    if(/Alignment\s+(\d+):/)
    {
        $ALN_ID = $1;
        $info{$ALN_ID} = $_;
        $strand{$ALN_ID} = "minus" if (/minus/);
        $strand{$ALN_ID} = "plus" if (/plus/);
    }
    else
    {
        chomp; 
        my @e=split;
        $start{$ALN_ID} = $e[1] if (!exists($start{$ALN_ID}));
        $end{$ALN_ID} = $e[1];
        $ortholog{$ALN_ID}{$e[0]} = $e[1];
        $ortholog{$ALN_ID}{$e[1]} = $e[0];
    }
}


my $this_gene;
for my $ALN_ID (sort {$a<=>$b} keys %start)
{
    $this_gene = $start{$ALN_ID};
    print $info{$ALN_ID};
    if($strand{$ALN_ID} eq "plus")
    {
        while( $next{$this_gene} ne $end{$ALN_ID})
        {
            print $this_gene, "   ";
            print "\t", $ortholog{$ALN_ID}{$this_gene} if (exists $ortholog{$ALN_ID}{$this_gene});
            print "\n";
            if(exists($next{$this_gene}))
            {
                ($this_gene) = ($next{$this_gene});
            }
            else
            {
                last;
            }
        }
    }
    elsif($strand{$ALN_ID} eq "minus")
    {
        while( $prev{$this_gene} ne $end{$ALN_ID})
        {
            print $this_gene, "   "; 
            print "\t", $ortholog{$ALN_ID}{$this_gene} if (exists $ortholog{$ALN_ID}{$this_gene});
            print "\n";
            if(defined $this_gene and exists($prev{$this_gene}))
            {
                ($this_gene) = ($prev{$this_gene}) ;
            }
            else
            {
                last;
            }
        }
    }
}
