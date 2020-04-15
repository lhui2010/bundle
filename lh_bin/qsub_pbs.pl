#!/usr/bin/perl -w


if(@ARGV<1)
{
	print STDERR  "Usage: $0 shell\n";
	print STDERR "I'll split shell into 20 partition and qsub to small que and wait for those jobs to end\n";
	exit;
}

$shell = $ARGV[0];

`perl /lustre/user/liuhui/bin/get_pbs_shell.pl $shell`;
`chmod 755 $shell.*.pbs`;
$list = `for i in $shell.*.pbs; do qsub -p 1023 -q long \$i; done`;

@list = split/\n/, $list;
for(@list)
{
	s/\..*//;
#	$hash{$_} = 1;
}

&wait_all_pbs;

`rm $shell.*.pbs`;


sub wait_all_pbs{
        while(1)
        {
		$mark = 1;
		$sys = `qstat -u liuhui|grep -v "C"`;
		for(@list)
		{
			if($sys=~/$_/)
			{
				 sleep 60;
				 $mark = 0;
				 last;
			}
		}
		if ($mark)
		{
			last;
		}
		
	}

}


