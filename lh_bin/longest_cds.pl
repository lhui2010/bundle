#!/usr/bin/perl -w

#gene:Os08g0254300

my $fname = $ARGV[0];
$/ = '>';
while(<>)
{
	my ($id, $seq) = split /\n/, $_, 2;

	next unless (defined $seq and $seq ne "");

    if($id=~/(^.*?) .* (gene:.*?) /)
    {
        $rna_id = $1;
        $gene_id = $2;
    }
    else
    {
        $id =~ s/\s.*//;
        #$id =~ s/\./_/g;
        $rna_id = $id;
        $gene_id = $id;
        $gene_id =~ s/T\d+$//;
        $gene_id =~ s/-t\d+$//;
        $gene_id =~ s/\.m\d+$//;
    }

	$seq =~ s/\n//g;
	$seq =~ s/\s//g;

	$len = length($seq);

    if(!exists $len_dict{$gene_id} or $len > $len_dict{$gene_id})
    {
        $len_dict{$gene_id} = $len;
        $gene_dict{$gene_id} = $rna_id;
    }
}

for my $k (sort keys %gene_dict)
{
    print $gene_dict{$k}, "\t", $k, "\n";
}


