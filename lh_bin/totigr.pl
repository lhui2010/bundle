
open IN, "RAP-MSU.txt" or die;

while(<IN>)
{
	@e=split;
	$rap_id = shift @e;
	$tigr_id = shift @e;
	$tigr_id =~ s/,//g if (defined $tigr_id);

	$converntor{$rap_id} = $tigr_id if (defined $tigr_id);
}

while(<>)
{
	chomp;
	print $converntor{$_},"\n" if (exists $converntor{$_});
}
	
