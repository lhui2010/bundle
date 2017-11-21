use strict;

$_=<>;
chomp;
my  @e=split;
my ($lid, $ltig, $lstart, $lend)  = ($e[0], $e[1], $e[2], $e[3]); #sorted
my $lprint = $_;
while(<>)
{
    chomp;
    my @e  = split;
    my ($id, $tig, $start, $end) = ($e[0], $e[1], $e[2], $e[3]);

    ($start, $end) = ($end, $start) if ($start > $end);

    if ($tig ne $ltig or $start > $lend)
    {
        print $lprint,"\n";
        ($lid, $ltig, $lstart, $lend)  = ($id, $tig, $start, $end);
        $lprint = $_;
    }
    elsif($start < $lend and $end >$lend)
    {
        $lend=$end;
        warn $.;
    }
}

print $_,"\n";
