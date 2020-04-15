#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage=<<INFO;
Description:
	2010-03-01
Author:
	Hezengquan, HSB
Usage: 
	perl $0 -f11 <file1.out> -f12 <file2.out> -f21 <file1.2seg> -f22 <file2.2seg> -o <junction file prefix> -i <insert length>
INFO


my ($file11,$file12,$file21,$file22,$head,$insertLength,$readLength,$help);
my ($keyDir,$keyReadID,$keyReadNum,$keyReadPos,$keyReadIDNoPE,$keyReadNumNoPE,$keyReadPosNoPE,$keySupportReadID,$keySupportReadNum)=('juncDir','readID','readNum','readPos', 'readIDNoPE', 'readNumNoPE', 'readPosNoPE', 'supportReadID','supportReadNum');
GetOptions(
#	"psl:s"=>\$psl,
	"f11:s"=>\$file11,
	"f12:s"=>\$file12,
	"f21:s"=>\$file21,
	"f22:s"=>\$file22,
	"o:s"=>\$head,
	"i:i"=>\$insertLength,
	"h|help"=>\$help
);

die "$usage" if(!$file11||!$file12||!$file21||!$file22||!$head||$help); 

my ($resultJUNC) = ($head.'.junc');

my %position;

my %oneSeg;
my %twoSeg;

my %myjunc;

my $startTime = time();
my $PESupport = 0;
ReadOneSeg($file11,'a',\%oneSeg);
ReadOneSeg($file12,'b',\%oneSeg);

ReadTwoSeg($file21,'a',\%twoSeg);
ReadTwoSeg($file22,'b',\%twoSeg);
my ($insertGapLB,$insertGapUB)=($insertLength*4/5-$readLength*2-10,$insertLength*6/5-$readLength*2+2000);

print STDOUT "\n";

open my $OUTJUNC,'>',$resultJUNC || die "Can't create $resultJUNC!";

DealTwoSeg($file21,'a');
DealTwoSeg($file22,'b');


CheckPairEnd();
      foreach my $i(sort keys %myjunc)
       {
               foreach my $j (sort {$a<=>$b} keys %{$myjunc{$i}})
               {
                       foreach my $k (sort {$a<=>$b} keys %{$myjunc{$i}{$j}})
                       {
			       print $OUTJUNC $i."\t".$j."\t".$k."\t".$myjunc{$i}{$j}{$k}{$keyDir}."\t";
			       print $OUTJUNC ($myjunc{$i}{$j}{$k}{$keyReadNum}+$myjunc{$i}{$j}{$k}{$keyReadNumNoPE}).":";
			       if(exists $myjunc{$i}{$j}{$k}{$keyReadNum})
			       {
				       $PESupport = $myjunc{$i}{$j}{$k}{$keyReadNum};
				       print $OUTJUNC $PESupport."|";
			       }
			       else
			       {
				       $PESupport = 0;
				       print $OUTJUNC "0|";
			       }
			       if(exists $myjunc{$i}{$j}{$k}{$keyReadNumNoPE})
			       {
				       print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadNumNoPE}."\t";
			       }
			       else
			       {
				       print $OUTJUNC "0\t";
			       }
			       #print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadNum}.":".$myjunc{$i}{$j}{$k}{$keyReadNumNoPE}."\t";
			       if(exists $myjunc{$i}{$j}{$k}{$keySupportReadNum})
			       {
				       $PESupport += $myjunc{$i}{$j}{$k}{$keySupportReadNum};
				       print $OUTJUNC $myjunc{$i}{$j}{$k}{$keySupportReadNum}."\t";
			       }
			       else
			       {
				       print $OUTJUNC "0\t";
			       }
			      	
			       if($PESupport>0)
			       {
				       print $OUTJUNC "PE\t";
			       }
			       else
			       {
				       print $OUTJUNC "SE\t";
			       }
			       if(exists $myjunc{$i}{$j}{$k}{$keyReadID})
			       {
				       print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadID}."|";
			       }
			       else
			       {
				       print $OUTJUNC "*|";
			       }
			       
			       if(exists $myjunc{$i}{$j}{$k}{$keyReadIDNoPE})
			       {
				       print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadIDNoPE}."\t";
			       }
			       else
			       {
				       print $OUTJUNC "*\t";
			       }
			       #print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadID}.":".$myjunc{$i}{$j}{$k}{$keyReadIDNoPE}."\t";

			       if(exists $myjunc{$i}{$j}{$k}{$keySupportReadID})
			       {
				       print $OUTJUNC $myjunc{$i}{$j}{$k}{$keySupportReadID}."\t";
			       }
			       else
			       {
				       print $OUTJUNC "*\t";
			       }

			       if(exists $myjunc{$i}{$j}{$k}{$keyReadPos})
			       {
				       print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadPos}."|";
			       }
			       else
			       {
				       print $OUTJUNC "*|";
			       }
			       
			       if(exists $myjunc{$i}{$j}{$k}{$keyReadPosNoPE})
			       {
				        print $OUTJUNC $myjunc{$i}{$j}{$k}{$keyReadPosNoPE}."\n";
			       }
			       else
			       {
				       print $OUTJUNC "*\n";
			       }
                       }
               }
      }

print STDOUT "Finish!\n";
my $time=time()-$startTime;
print STDOUT "Elapsed time:\t".$time."s\n";


sub DealTwoSeg{
	my ($file,$tag)=@_;
	my $untag = ($tag eq 'a')? 'b':'a';
	my $readTag = ($tag eq 'a') ? '1':'2';

	open STDIN,'<',$file || die "Can't open $file!";
	while (<STDIN>) {
		next if(/^\s*$/);
		chomp;
		my @temp = split(/\t/,$_);
		next if($temp[2] ne "1,1");
		my @chrID = split(/,/,$temp[6]);
		next if($chrID[0] ne $chrID[1]);
		
		my $strand = $temp[4];
		if($twoSeg{$temp[0]}{$tag}{'times'}>=2 && $twoSeg{$temp[0]}{$tag}{'strand'} ne $twoSeg{$temp[0]}{$tag}{'strand2'})
		{
			if(exists($oneSeg{$temp[0]}{$untag}) && $oneSeg{$temp[0]}{$untag} ne ''){
				
				if(substr($oneSeg{$temp[0]}{$untag},0,1) eq substr($strand,0,1)) # HSB add[0]
				{ 
					next;
				}
			}
			elsif($twoSeg{$temp[0]}{$untag}{'times'}>=2 && $twoSeg{$temp[0]}{$untag}{'strand'} ne $twoSeg{$temp[0]}{$untag}{'strand2'})
			{
				next;
			}
			elsif(exists($twoSeg{$temp[0]}{$untag}) && $twoSeg{$temp[0]}{$untag}{'strand'} ne ''){
				if($twoSeg{$temp[0]}{$untag}{'strand'} eq $strand){
					next;
				}
			}
			else
			{
				next;
			}
		}

		my @len = split(/,/,$temp[3]);
		my @start = split(/,/,$temp[7]);

		my ($JuncPos1,$JuncPos2);
		$JuncPos1 = $start[0]+$len[0]-1;
		$JuncPos2 = $start[1];

		next if($JuncPos2-$JuncPos1-1>=500000);
		next if($JuncPos2-$JuncPos1-1<=50);
		
		my $direction = $temp[5];
		my $PEExist = &CheckPairEndExist($temp[0],$temp[4],$chrID[0],$temp[7],$temp[3],$tag,$direction); #readID, strand, chrID, startPos[2],length[2],tag(a or b)
		next if($PEExist == -1);
		$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyDir}=$direction;
		if($PEExist == 1)
		{
			#my %junc=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2};
			if(exists $myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadID})
			{
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadID} = $myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadID}.",".$readTag."-".$temp[0];
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPos} = $myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPos}.",".$start[0];
			}
			else
			{
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadID} = $readTag."-".$temp[0];
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPos} = $start[0];			
			}
			#$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadID}=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadID}.$readTag."-".$temp[0].",";
			#$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPos}=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPos}.$start[0].",";
			$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadNum}++;
		}
		elsif($PEExist==0) #$PEExist==-1, don't treat it as a junction
		{
			if(exists $myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadIDNoPE})
			{
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadIDNoPE}=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadIDNoPE}.",".$readTag."-".$temp[0];
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPosNoPE}=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPosNoPE}.",".$start[0];
			}
			else
			{
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadIDNoPE}=$readTag."-".$temp[0];
				$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPosNoPE}=$start[0];
			}
			#$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadIDNoPE}=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadIDNoPE}.$readTag."-".$temp[0].",";
			#$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPosNoPE}=$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadPosNoPE}.$start[0].",";
			$myjunc{$chrID[0]}{$JuncPos1}{$JuncPos2}{$keyReadNumNoPE}++;
		}
	}
	close STDIN;
}
sub CheckPairEndExist
{
	my ($readID, $strand, $chrID, $positions, $lengths, $tag, $direction) = @_;
	my @pos = split (/,/,$positions);
	my @len = split (/,/,$lengths);
	my $untag = ($tag eq 'a')? 'b':'a';
	my ($left,$right,$gap);
	if($twoSeg{$readID}{$tag}{'times'}>1)
	{
		return 1;
	}
	#PE is one segment
	if(exists $oneSeg{$readID}{$untag}&& $oneSeg{$readID}{$untag} ne '')
	{
		my @oneSegInfo = split (/,/,$oneSeg{$readID}{$untag});
		my ($strandOne, $chrIDOne, $posOne, $lengthOne) = ($oneSegInfo[0],$oneSegInfo[1],$oneSegInfo[2],$oneSegInfo[3]);
		if($strandOne eq $strand)
		{
			return -1;
		}
		if($strandOne ne $strand && $chrIDOne eq $chrID)
		{
			if($posOne < $pos[0])
			{
				$left = $posOne + $lengthOne -1;
				$right = $pos[0];
			}
			else
			{
				$left = $pos[1] + $len[1] - 1;
				$right = $posOne;
			}
			$gap = $right - $left - 1;
		if($gap>=$insertGapLB && $gap<=$insertGapUB)
			{
				return 1;
			}
		}
	}
	elsif(exists $twoSeg{$readID}{$untag} && $twoSeg{$readID}{$untag}{'strand'} ne "")
	{
		if($twoSeg{$readID}{$untag}{'times'}==1)
		{
			if($twoSeg{$readID}{$untag}{'strand'} ne $strand)
			{
				my @twoSegResult = split (/,/,$twoSeg{$readID}{$untag}{'alignment'});
				my $chrIDTwo = $twoSegResult[0];
				my ($posOne,$posTwo,$lenOne,$lenTwo,$directionTwo) = ($twoSegResult[1],$twoSegResult[2],$twoSegResult[3],$twoSegResult[4],$twoSegResult[5]);
				
				if($chrID eq $chrIDTwo && $direction eq $directionTwo)
				{
					if($posOne < $pos[0])
					{
						$left = $posTwo + $lenTwo -1;
						$right = $pos[0];
					}
					else
					{
						$left = $pos[1] + $len[1] -1;
						$right = $posOne;
					}
					$gap = $right - $left - 1;
	                        if($gap>=$insertGapLB && $gap<=$insertGapUB)
					{
						return 1;
					}
				}
				else
				{
					return -1;
				}
			}
		}
		elsif($twoSeg{$readID}{$untag}{'times'}>1)
		{
			if($twoSeg{$readID}{$untag}{'strand'} ne $strand)
			{
				my @twoSegResult = split (/,/,$twoSeg{$readID}{$untag}{'alignment'});
				my $chrIDTwo = $twoSegResult[0];
				my ($posOne,$posTwo,$lenOne,$lenTwo,$directionTwo) = ($twoSegResult[1],$twoSegResult[2],$twoSegResult[3],$twoSegResult[4],$twoSegResult[5]);
				if($chrID eq $chrIDTwo && $direction eq $directionTwo)
				{
					if($posOne < $pos[0])
					{
						$left = $posTwo + $lenTwo -1;
						$right = $pos[0];
					}
					else
					{
						$left = $pos[1] + $len[1] -1;
						$right = $posOne;
					}
					$gap = $right - $left - 1;
	                        if($gap>=$insertGapLB && $gap<=$insertGapUB)
					{
						return 1;
					}
				}
				else
				{
					return -1;
				}
					
			}
			elsif($twoSeg{$readID}{$untag}{'strand2'} ne $strand)
			{
				my @twoSegResult = split (/,/,$twoSeg{$readID}{$untag}{'alignment2'});
				my $chrIDTwo = $twoSegResult[0];
				my ($posOne,$posTwo,$lenOne,$lenTwo,$directionTwo) = ($twoSegResult[1],$twoSegResult[2],$twoSegResult[3],$twoSegResult[4],$twoSegResult[5]);
				if($chrID eq $chrIDTwo && $direction eq $directionTwo)
				{
					if($posOne < $pos[0])
					{
						$left = $posTwo + $lenTwo -1;
						$right = $pos[0];
					}
					else
					{
						$left = $pos[1] + $len[1] -1;
						$right = $posOne;
					}
					$gap = $right - $left - 1;
	                        if($gap>=$insertGapLB && $gap<=$insertGapUB)
					{
						return 1;
					}
				}
				else
				{
					return -1;
				}
			}
		}
	}
	return 0;
}
sub CheckPairEnd
{
	foreach my $readID (keys %oneSeg)
	{
		if(exists $oneSeg{$readID}{'a'})
		{
			if(exists $oneSeg{$readID}{'b'})
			{
				&CheckOneOnePairEnd($readID,$oneSeg{$readID}{'a'},$oneSeg{$readID}{'b'});
			}
			elsif(exists $twoSeg{$readID}{'b'} && $twoSeg{$readID}{'b'}{'strand'} ne '')
			{
				my $twoSegRead = $twoSeg{$readID}{'b'}{'strand'}.",".$twoSeg{$readID}{'b'}{'alignment'}; 
				&CheckOneTwoPairEnd($readID,$oneSeg{$readID}{'a'},$twoSegRead);				
				if($twoSeg{$readID}{'b'}{'times'}>=2)
				{
					$twoSegRead = $twoSeg{$readID}{'b'}{'strand2'}.",".$twoSeg{$readID}{'b'}{'alignment2'};
					&CheckOneTwoPairEnd($readID,$oneSeg{$readID}{'a'},$twoSegRead);
				}
			}
		}
	}
	
	foreach my $readID (keys %twoSeg)
	{
		if(exists $twoSeg{$readID}{'a'} && $twoSeg{$readID}{'a'}{'strand'} ne '')
		{
			my $twoSegReadOne = $twoSeg{$readID}{'a'}{'strand'}.",".$twoSeg{$readID}{'a'}{'alignment'};
			#print "'a' 2-seg,ReadID:".$readID."\n";
			if(exists $oneSeg{$readID}{'b'})
			{
				&CheckOneTwoPairEnd($readID,$oneSeg{$readID}{'b'},$twoSegReadOne);
				#print "ReadID:21:".$readID."\n";
				if($twoSeg{$readID}{'a'}{'times'}>=2)
				{
					$twoSegReadOne = $twoSeg{$readID}{'a'}{'strand2'}.",".$twoSeg{$readID}{'a'}{'alignment2'};
					&CheckOneTwoPairEnd($readID,$oneSeg{$readID}{'b'},$twoSegReadOne);
				}
			}
			elsif(exists $twoSeg{$readID}{'b'} && $twoSeg{$readID}{'b'}{'strand'} ne '')
			{
				my $twoSegReadTwo = $twoSeg{$readID}{'b'}{'strand'}.",".$twoSeg{$readID}{'b'}{'alignment'};
				&CheckTwoTwoPairEnd($readID, $twoSegReadOne, $twoSegReadTwo);
				if($twoSeg{$readID}{'a'}{'times'}>=2)
				{
				        $twoSegReadOne = $twoSeg{$readID}{'a'}{'strand2'}.",".$twoSeg{$readID}{'a'}{'alignment2'};
				        &CheckTwoTwoPairEnd($readID, $twoSegReadOne, $twoSegReadTwo);
				}
				elsif($twoSeg{$readID}{'b'}{'times'}>=2)
				{
					$twoSegReadTwo = $twoSeg{$readID}{'b'}{'strand2'}.",".$twoSeg{$readID}{'b'}{'alignment2'};
					&CheckTwoTwoPairEnd($readID, $twoSegReadOne, $twoSegReadTwo);
				}
			}	
		}
	}
}

sub CheckOneOnePairEnd
{
	my ($readID,$valueOne,$valueTwo) = @_;
	my @infoOne = split(/,/,$valueOne);
	my @infoTwo = split(/,/,$valueTwo);
	my ($strandOne,$chrIDOne,$posOne,$lengthOne) = ($infoOne[0],$infoOne[1],$infoOne[2],$infoOne[3]);
	my ($strandTwo,$chrIDTwo,$posTwo,$lengthTwo) = ($infoTwo[0],$infoTwo[1],$infoTwo[2],$infoTwo[3]);
	if($strandOne ne $strandTwo && $chrIDOne eq $chrIDTwo)
	{
		my ($left,$right);
		if($posOne<$posTwo)
		{
			$left = $posOne+$lengthOne-1;
			$right = $posTwo;
		}
		else
		{
			$left = $posTwo+$lengthTwo-1;
			$right = $posOne;
		}
		&CheckSupport($readID,$chrIDOne,$left,$right);
	}
}

sub CheckOneTwoPairEnd
{
	my ($readID,$oneSegResult,$twoSegResult) = @_;
	my @oneSegInfo = split(/,/,$oneSegResult);
	my @twoSegInfo = split(/,/,$twoSegResult);
	my ($strandOne, $chrIDOne, $posOne, $lengthOne) = ($oneSegInfo[0],$oneSegInfo[1],$oneSegInfo[2],$oneSegInfo[3]);
	my ($strandTwo, $chrIDTwo, $twoSegPosOne,$twoSegPosTwo,$twoSegLengthOne,$twoSegLengthTwo) = ($twoSegInfo[0],$twoSegInfo[1],$twoSegInfo[2],$twoSegInfo[3],$twoSegInfo[4],$twoSegInfo[5]);
	if($strandOne ne $strandTwo && $chrIDOne eq $chrIDTwo)
	{
		my ($left,$right);
		if($posOne<$twoSegPosOne)
		{
			$left = $posOne+$lengthOne-1;
			$right = $twoSegPosOne;
			&CheckSupport($readID,$chrIDOne,$left,$right);
		}
		elsif($posOne>$twoSegPosTwo)
		{
			$left = $twoSegPosTwo + $twoSegLengthTwo - 1;
			$right = $posOne;
			&CheckSupport($readID,$chrIDOne,$left,$right);
		}
	}
}

sub CheckTwoTwoPairEnd
{
	my ($readID, $twoSegResultOne,$twoSegResultTwo) = @_;
	my @twoSegInfoOne = split(/,/,$twoSegResultOne);
	my @twoSegInfoTwo = split(/,/,$twoSegResultTwo);
	my ($strandOne,$chrIDOne,$twoSegPosOne,$twoSegPosTwo,$twoSegLengthOne,$twoSegLengthTwo) = ($twoSegInfoOne[0],$twoSegInfoOne[1],$twoSegInfoOne[2],$twoSegInfoOne[3],$twoSegInfoOne[4],$twoSegInfoOne[5]);
	my ($strandTwo,$chrIDTwo,$twoSegPosThree,$twoSegPosFour,$twoSegLengthThree,$twoSegLengthFour) = ($twoSegInfoTwo[0],$twoSegInfoTwo[1],$twoSegInfoTwo[2],$twoSegInfoTwo[3],$twoSegInfoTwo[4],$twoSegInfoTwo[5]);
	if($strandOne ne $strandTwo && $chrIDOne eq $chrIDTwo)
	{
		my ($left, $right);
		if($twoSegPosTwo<$twoSegPosThree)
		{
			$left = $twoSegPosTwo + $twoSegLengthTwo -1;
			$right = $twoSegPosThree;
			&CheckSupport($readID,$chrIDOne,$left,$right);
		}
		elsif($twoSegPosFour<$twoSegPosOne)
		{
			$left = $twoSegPosFour + $twoSegLengthFour -1;
			$right = $twoSegPosOne;
			&CheckSupport($readID,$chrIDOne,$left,$right);
			
		}
	}
}

#argv: $readID,$chrID, $left,$right, check which junction it supports
sub CheckSupport
{
	my ($readID,$chrID,$left,$right) = @_;
	my $i;
	#my ($insertGapLB,$insertGapUB)=($insertLength*4/5-$readLength*2,$insertLength*6/5-$readLength*2);
	my $tempInsertGapUB = $insertLength*6/5-$readLength*2 + 100;
	if(exists $myjunc{$chrID})
	{
		for($i=0;$i<=$tempInsertGapUB; $i++)
		{
			if(exists $myjunc{$chrID}{$left+$i})
			{
				#print $chrID."\n";
				#print $left+$i." exists\n";
				foreach my $rightJunc (keys %{$myjunc{$chrID}{$left+$i}})
				{
					if($rightJunc <= $right && $i+$right-$rightJunc<=$tempInsertGapUB)
					{
						#$myjunc{$chrID}{$left+$i}{$rightJunc}{"3-".$readID}=1;
						if(exists $myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadID})
						{
							$myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadID}=$myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadID}.",".$readID;
						}
						else
						{
							$myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadID}=$readID;
						}
						#$myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadID}=$myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadID}.$readID.",";
						$myjunc{$chrID}{$left+$i}{$rightJunc}{$keySupportReadNum}++;
					}
					
				}
			}
		}
	}
}

# read 1 segment file
sub ReadOneSeg{
	my ($file,$tag,$hash) = @_;
	open STDIN,'<',$file || die "Can't open $file!";
	while (<STDIN>) {
		next if(/^\s*$/);
		chomp;
		my @lines = split /\t/,$_;
		next if($lines[4] ne "1");

		my $chrID = $lines[7];
		my $pos = $lines[8];
		my $length = $lines[5];
		my $strand = $lines[6];
		$hash->{$lines[0]}{$tag} = $strand.",".$chrID.",".$pos.",".$length;
	}
	close STDIN;
	print STDOUT "finish reading $file!\n";
}

# read 2 segment file
sub ReadTwoSeg{
	my ($file,$tag,$hash) = @_;
	my @lines;
	open STDIN,'<',$file || die "Can't open $file!";
	while (<STDIN>) {
		next if(/^\s*$/);
		chomp;
		@lines = split /\t/,$_;
		next if($lines[2] ne "1,1");

		my @chrID = split(/,/,$lines[6]);
		next if($chrID[0] ne $chrID[1]);
		
		if(exists $hash->{$lines[0]}{$tag}{'times'} &&$hash->{$lines[0]}{$tag}{'times'}==2 )
		{
			$hash->{$lines[0]}{$tag}{'strand3'} = $lines[4];
			$hash->{$lines[0]}{$tag}{'alignment3'} = $chrID[0].",".$lines[7].",".$lines[3].",".$lines[5];
		}
		elsif(exists $hash->{$lines[0]}{$tag}{'times'})
		{
			$hash->{$lines[0]}{$tag}{'strand2'} = $lines[4];
			$hash->{$lines[0]}{$tag}{'alignment2'} = $chrID[0].",".$lines[7].",".$lines[3].",".$lines[5];
			                                        # chrID,postions,lengths,direction
		}
		else
		{
			$hash->{$lines[0]}{$tag}{'strand'} = $lines[4];
			$hash->{$lines[0]}{$tag}{'alignment'} = $chrID[0].",".$lines[7].",".$lines[3].",".$lines[5];
		}
		$hash->{$lines[0]}{$tag}{'times'} ++;			
	}
	my @len = split(/,/,$lines[3]);
	$readLength = $len[0] + $len[1];
	close STDIN;
	print STDOUT "finish reading $file!\n";
}

