#B73V4_ctg98 gramene exon    36196   36270   .   -   .   Parent=transcript:Zm00001d000230_T002
#

my %chr_id;
while(<>)
{
    my ($seq_id, undef, undef, $start, $end, undef, $strand, undef, $transcript_id)= split;
    ($end, $start) = ($start, $end) if($strand eq "-");

    #warn $seq_id; exit;

    $transcript_id=~s/Parent=transcript://;

	$strand{$transcript_id} = $strand;
    $count{$transcript_id}++; 
    $chr_id{$transcript_id} = $seq_id;
    if($count{$transcript_id} == 1)
    {   
        $hash{$transcript_id}="Esngl\t$start\t$end\t$transcript_id\n"
    }
    else
    {   
        $hash{$transcript_id}=~s/Esngl/Einit/;
        $hash{$transcript_id}.="Exon\t$start\t$end\t$transcript_id\n"
    }
}

#for my $k (sort keys %hash)
#for my $k (sort {$chr_id{$a} cmp $chr_id{$a}} keys %hash)
my %print;
for my $k (sort {$chr_id{$a} cmp $chr_id{$b}} keys %chr_id)
{
    if(!exists($print{$chr_id{$k}}))
    {   
        print ">$chr_id{$k}\n";
        $print{$chr_id{$k}} = 1;
    }
    (my $result = $hash{$k}) =~ s{(.*)Exon}{$1Eterm}xms;
	
	if($strand{$k} eq "-")
	{
		my @e=split/\n/, $result;
		@re=reverse @e;
		$re[0]=~s/Eterm/Einit/;
		$re[-1]=~s/Einit/Eterm/;
		$result = join ("\n", @re);
		$result.="\n";
	}
	
	
    print "$result";
}
