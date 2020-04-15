#!/usr/bin/perl -w

print ".pl fa (chromosome) start end >result
seq[0]=>
"if(@ARGV==0);

$seq=">";


open IN, $ARGV[0] or die;

shift @ARGV;

$qry_chr = "";
if(@ARGV == 3)
{
	$qry_chr = shift @ARGV;
}

$start = shift @ARGV;
$end = shift @ARGV;



my $mark = 0;


my $rc_mark = 0;

if($end <$start)
{
	$rc_mark =1; 
	($end,$start) = ($start, $end);
}

while(<IN>)
{
	chomp;
	if(/>/)
	{
		if($qry_chr eq "")
		{
			$description=$_;
			$seq = ">";
			$mark = 1;
		}
		else
		{
			my @e=split;
			$chr_id = substr ($e[0], 1);		
			if($qry_chr eq $chr_id)
			{
				$description=$_;
				$seq = ">";
				$mark = 1;
			}
		}
	}
	elsif($mark==1)
	{
		$seq.=$_;
	}
}

#for $pos (@poslist)
#{
#	my ($start, $end) = @$pos;
	$seqout=substr($seq, $start, $end-$start+1);
	if($rc_mark)
	{
		&rc(\$seqout);
	}

	#print $description, "\t$start-$end\n";
	print $description, ":$start-$end\n";
	&formated_print($seqout);
#}

sub rc
{
    my $seq_p=shift;
    if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
        $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
        $$seq_p=reverse($$seq_p);
    }
}

sub formated_print
{
    my (@e1, $no);
    @e1=split//, $_[0];
    for($no=0; $no<@e1; $no++)
    {
        print $e1[$no];
        print "\n" if(($no+1)%70 == 0);
    }
    print "\n" unless ($no%70 == 0);
}
	
