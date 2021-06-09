#!/usr/bin/perl -w

my $column = shift;

#$_=<>;
my $count = 0;
my $sum = 0;
while(<>)
{
#	/size([0-9]*)/;
	my @e=split;
	if($e[$column]=~/^\d+$/)
    {
        $sum+=$e[$column];
        $count ++;
    }
}
print $sum/$count, "\n";
