while(<>)
{
	chomp;
	my ($qname, $sname, $identity, $match_len, $mis_len, undef, undef, undef, undef, undef, $evalue, $score) = split;

	$hash{$qname}{$sname} += $match_len;
	$print{$qname}{$sname} .= $_;
}

for my $k (sort keys %hash)
{
	for my $kk( sort {$hash{$k}{$b} <=> $hash{$k}{$a}} keys %{$hash{$k}})
	{
		print "$k\t$kk\t$hash{$k}{$kk}\t $print{$k}{$kk}\n";
		last;
	}
}
