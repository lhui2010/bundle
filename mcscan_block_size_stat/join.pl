
$_=<>;
chomp;
my ($lid, $ltig, $lstart, $lend)  = split; #sorted
$lprint = $_;
while(<>)
{
    chomp;
    my ($id, $tig, $start, $end)  = split;

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
    }
}

print $_,"\n";
