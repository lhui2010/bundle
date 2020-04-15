#!/usr/bin/perl -w

print "
this software is used to overcome shell's limitation on argument's size
usage:
.pl xx yy >result
where the filename is xx(*)yy
if xx is -, .pl will regard files' name as (*)yy
" if (@ARGV<1);
$prefix = $ARGV[0];
$prefix = "" if ($prefix eq "-");

@fnames=glob($prefix."*".$ARGV[1]);

for $fname(@fnames)
{
	open IN, $fname or die;
	while(<IN>)
	{
		print $_;
	}
}
