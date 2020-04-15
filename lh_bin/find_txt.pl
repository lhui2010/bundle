#!/usr/bin/perl -w
$key = shift @ARGV;
@fnames = @ARGV;

for $fname (@fnames)
{
	open IN, $fname or die;
	my $i = 0;
	while(<IN>)
	{
		$i++;
		if(/$key/)
		{
			print $fname,": $i\n";
			print $_;
		}
	}
}
