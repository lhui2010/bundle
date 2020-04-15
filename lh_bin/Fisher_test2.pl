#!/usr/bin/perl -w

$_=<>;
while(<>)
{
	chomp;
$line=$_;
@data=split(/\s+/,$line);
$id = shift @data;
$p_value=&fisher(@data);
print $id, "\t", $p_value,"\n";
}

sub fisher
{
	my ($cs,$us,$cn,$un)=@_;
	print "$cs\t$us\t$cn\t$un\t";
	my @matrix; my $N; my @R; my @C; 
	my $temp_z=0.0; my $temp_m=0.0; my $i; my $j; my $l;
	my $Pc=0.0; my $Pt=0.0; my $Ps=0.0;
	
	$matrix[0][0]=$cs;    $matrix[0][1]=$cn;
	$matrix[1][0]=$us;    $matrix[1][1]=$un;
	$R[0]=$cs+$cn; $R[1]=$us+$un;
	$C[0]=$cs+$us; $C[1]=$cn+$un;
	$N=$R[0]+$R[1];
	
	for($i=0;$i<2;$i++){
		$temp_z+=&factorial($R[$i],0);
		$temp_z+=&factorial($C[$i],0);
		for($j=0;$j<2;$j++){
			$temp_m+=&factorial($matrix[$i][$j],0);
		}
	}
	$temp_m+=&factorial($N,0);
	$Pc=$temp_z-$temp_m;
	$Pc=exp($Pc);
	
	for($i=0;$i<=$R[0];$i++){
		$matrix_t[0][0]=$i;       $matrix_t[0][1]=$R[0]-$i;
		$matrix_t[1][0]=$C[0]-$i; $matrix_t[1][1]=$R[1]-$C[0]+$i;
		$temp_z=0.0;              $temp_m=0.0;
		next if ($matrix_t[0][0]<0 || $matrix_t[0][1]<0 || $matrix_t[1][0]<0 || $matrix_t[1][1]<0);
		for($l=0;$l<2;$l++){
			$temp_z+=&factorial($R[$l],0);
			$temp_z+=&factorial($C[$l],0);
			for($j=0;$j<2;$j++){
				$temp_m+=&factorial($matrix_t[$l][$j],0);
			}
		}
		$temp_m+=&factorial($N,0);
		$Pt=$temp_z-$temp_m;
		$Pt=exp($Pt);
		if($Pt<=$Pc){
			$Ps+=$Pt;
		}
	}
	return ($Ps);
}

sub factorial
{
	my ($n,$num)=@_;
	my $i; my $temp=1; my $seta=0.0;
	
	if($n<0){
		print "ERROR in factorial!\n";
		exit;
	}
	$n = $n + 1;
		$x = 0;
                $x += 0.1659470187408462e-06/($n+7);
                $x += 0.9934937113930748e-05/($n+6);
                $x -= 0.1385710331296526    /($n+5);
                $x += 12.50734324009056     /($n+4);
                $x -= 176.6150291498386     /($n+3);
                $x += 771.3234287757674     /($n+2);
                $x -= 1259.139216722289     /($n+1);
                $x += 676.5203681218835     /($n);
                $x += 0.9999999999995183;
                $temp = log($x) - 5.58106146679532777 - $n +($n - 0.5)*log($n + 6.5);
	return ($temp);
}

