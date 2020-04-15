#!/usr/bin/perl -w

#add deredundant func
#12:34 2010/8/1

#used to count aligned reads number from accepted hits .sam file
#9:55 2010/8/1

$file=$ARGV[0];

$current_line=<>;
$current_line=<>;
$current_line=<>;
@e=split /\t/, $current_line;
$current_name{$e[0]}=1;

while(<>)
{
    chomp;
	@e=split /\t/, $_;
	$current_name{$e[0]}=1;
}

@namelist=sort keys %current_name;

$size=@namelist;

print "$file\t$size\n";
