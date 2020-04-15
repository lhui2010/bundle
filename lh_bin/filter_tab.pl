$_=<>; print ;
while(<>)
{
	$count =0;
	chomp;
	@e=split;
	my $id = shift @e;shift @e;
	$sp_num = $#e;
	pop @e;

	for my $i(@e)
	{
		if ($i > 0)
		{
			$count ++;
		}
	}
	if($sp_num - $count <=2)
	{
		$_ = join "\t", ($id, $id, @e); 
		print $_,"\n";
	}
}
