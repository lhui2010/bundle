#!/usr/bin/perl -w
$line = `ps x`, ;
@lines = split /\n/, $line;
$i = 0;

for (@lines)
{
	next if (/sleep.pl/);
	$i++ if (/blat/ or /perl/ or /gmhmme3/);
}

while($i > 0)
{
	sleep(10);
	
	$i=0;
	$line = `ps x`;
	@lines = split /\n/, $line;
	
	for (@lines)
	{
		next if (/sleep.pl/);
	   $i++ if (/blat/ or /perl/ or /gmhmme3/);
	}
}
	
