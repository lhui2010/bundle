#!/usr/bin/perl -w


$thread = 60;
$dir = `pwd`;
$count = 1;
$fname = $ARGV[-1];
$total = `wc -l $fname`;
($total) = (split /\s+/, $total)[0];

$max_line_per_file = int($total/$thread);

open OUT,">$fname.$count.pbs";

print OUT "cd $dir\n";
$count_line = 1;
while(<>)
{
	$count_line++;
#	print OUT "cd $dir\n";
	print OUT $_, "\n";
	if($count_line >$max_line_per_file)
	{
		$count_line = 1;
		if ($.>=$total)
		{
			exit;
		}
		close OUT;
		$count++;
		open OUT,">$fname.$count.pbs";
		print OUT "cd $dir\n";
	}
}


