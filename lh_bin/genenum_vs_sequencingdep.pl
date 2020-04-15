#!/usr/bin/perl -w

#only look for anovel miRNA, ignoring ncRNAxxx
#2010/8/11

#usage: thisprogram.pl *.in_Genome
#result will be written to origo file in this format:
#cellname uniqreadscount(1..10) real count

#fixed som bugs
#2010/8/6



#will add a new func, total uniq count
#18:30 2010/8/5

@filenames=@ARGV;

open Fout, ">origo" || die "can't open output \n";

for $filename(@filenames)
{
	open Fin, $filename || die "can't open input \n";

	@splitfilename=split/\./, $filename;
	$cellname=$splitfilename[0];

	print Fout "$cellname\t";

	undef @purenamelist;
	undef @namelist;

	while(<Fin>)
	{
		@elem=split;
		$name=$elem[0];
		$count=$elem[1];

		next if ( substr($name, 0, 5) eq "ncRNA");# excluding novel small RNA, which may be introns.		

		push(@purenamelist, $name);#no extra redundance

		for(1..$count)
		{
			push @namelist, $name;#full redundance
		}
	}

	for($scale=1; $scale<=10; $scale++)
	{
		undef @rndlist;
		undef @uniq;

		$size=int (@namelist*$scale/10);
		for($i=0; $i<$size; $i++)
		{
			push @rndlist, $namelist[int(rand(@namelist))];
		}
		@sorted = sort @rndlist;

		$last_elem = 0;

		for $sort_elem(@sorted)
		{
			push @uniq, $sort_elem if($sort_elem ne $last_elem);
			$last_elem = $sort_elem;
		}
		$size=@uniq;
		print Fout "$size\t";
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
