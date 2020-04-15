#!/usr/bin/perl -w

#will
#transform


##left column is unique, right column contains duplicated ids
#Os01t0217300    IPR000169
#Os01t0347600    IPR000169
#Os01t0612700    IPR000169

#to

#IPR000169 Os01t0217300 Os01t0347600 Os01t0612700


my %id_collection;
my %id_gene_count;
while(<>)
{
	my ($gene, @id) = split /\s+/, $_;

	next if (@id < 1);# some genes do not have a go id

	for my $id (@id)
	{
		$id_collection{$id} .= $gene."\t";
		$id_gene_count{$id} ++;
	}
}

for my $key (sort { $id_gene_count{$b} <=> $id_gene_count{$a}} keys %id_gene_count )
{
	print "$key\t$id_gene_count{$key}\t$id_collection{$key}\n";
}
