
$_=<>;
chomp;
my ($lid, $ltig, $lstart, $lend)  = split; #sorted
$lprint = $_;
while(<>)
{
    chomp;
    my ($id, $tig, $start, $end)  = split;

    if ($tig ne $ltig or $start > $lend)
    {
        print $lprint,"\n";
        ($lid, $ltig, $lstart, $lend)  = ($id, $tig, $start, $end);
        $lprint = $_;
    }

    if($start < $end and $end >$lend)
    {
        $lend=$end;
    }
}

print $_,"\n";
