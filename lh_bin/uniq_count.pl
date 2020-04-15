#!/usr/bin/perl -w

#usage£º ¡£¡£¡£pl *.in_Genome
#result will be written to origo file in this format:
#cellname uniqreadscount(1..10) real count

#fixed som bugs
#2010/8/6



#will add a new func, total uniq count
#18:30 2010/8/5

@filenames=@ARGV;

open Fout, ">uniq" || die "can't open output \n";

for $filename(@filenames)
{
	open Fin, $filename || die "can't open input \n";

	@splitfilename=split/\./, $filename;
	$cellname=$splitfilename[0];

	print Fout "$cellname\t";

	undef @purenamelist;

	while(<Fin>)
	{
		@elem=split;
		$name=$elem[0];
		push(@purenamelist, $name);#no extra redundance

	}


	#calculate real unique counts
	undef @uniq;
	@sorted = sort @purenamelist;
	$last_elem = 0;
	for $sort_elem(@sorted)
	{
		push @uniq, $sort_elem if($sort_elem ne $last_elem);
		$last_elem = $sort_elem;
	}
	$size=@uniq;
	print Fout "$size\n";
}

close Fout;
