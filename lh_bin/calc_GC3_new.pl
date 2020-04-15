#!/usr/bin/perl -w

#modified to calc GC3

#add the radical change count

#add gap count and length count func

#tested on one gene and two different genes, success
#getsnp.pl will still read R3, the ignorance will take place in this section
#change: ignore R3 for it's variable.success

#usage:getsnp45.pl genotype >result
@codon_name=qw/AAA AAT AAC AAG ATA ATT ATC ATG ACA ACT ACC ACG AGA AGT AGC AGG TAA TAT TAC TAG TTA TTT TTC TTG TCA TCT TCC TCG TGA TGT TGC TGG CAA CAT CAC CAG CTA CTT CTC CTG CCA CCT CCC CCG CGA CGT CGC CGG GAA GAT GAC GAG GTA GTT GTC GTG GCA GCT GCC GCG GGA GGT GGC GGG/;

#total family
my (%Counter_Total, %RSCU_Total, $Counter_GC3_Total, $Counter_Codon_Total);

my %number_codon;
for ($i = 0; $i <64; $i++)
{
	$number_codon{$codon_name[$i]} = $i;
}

push @{$aminoacid{"I"}}, qw/ATT ATC ATA/;
push @{$aminoacid{"L"}}, qw/CTT CTC CTA CTG TTA TTG/;
push @{$aminoacid{"V"}}, qw/GTT GTC GTA GTG/;
push @{$aminoacid{"F"}}, qw/TTT TTC/;
push @{$aminoacid{"M"}}, qw/ATG/;
push @{$aminoacid{"C"}}, qw/TGT TGC/;
push @{$aminoacid{"A"}}, qw/GCT GCC GCA GCG /;
push @{$aminoacid{"G"}}, qw/GGT GGC GGA GGG /;
push @{$aminoacid{"P"}}, qw/CCT CCC CCA CCG/;
push @{$aminoacid{"T"}}, qw/ACT ACC ACA ACG/;
push @{$aminoacid{"S"}}, qw/TCT TCC TCA TCG AGT AGC/;
push @{$aminoacid{"Y"}}, qw/TAT TAC/;
push @{$aminoacid{"W"}}, qw/TGG/;
push @{$aminoacid{"Q"}}, qw/CAA CAG/;
push @{$aminoacid{"N"}}, qw/AAT AAC/;
push @{$aminoacid{"H"}}, qw/CAT CAC/;
push @{$aminoacid{"E"}}, qw/GAA GAG/;
push @{$aminoacid{"D"}}, qw/GAT GAC/;
push @{$aminoacid{"K"}}, qw/AAA AAG/;
push @{$aminoacid{"R"}}, qw/CGT CGC CGA CGG AGA AGG/;
push @{$aminoacid{"*"}}, qw/TAA TAG TGA /;

#eg:input ATT, return I
sub nuc2pep
{
	return "I" if($_[0] eq "ATT" or $_[0] eq "ATC" or $_[0] eq "ATA");
	return "L" if($_[0] eq "CTT" or $_[0] eq "CTC" or $_[0] eq "CTA" or $_[0] eq "CTG" or $_[0] eq "TTA" or $_[0] eq "TTG");
	return "V" if($_[0] eq "GTT" or $_[0] eq "GTC" or $_[0] eq "GTA" or $_[0] eq "GTG");
	return "F" if($_[0] eq "TTT" or $_[0] eq "TTC");
	return "M" if($_[0] eq "ATG");
	return "C" if($_[0] eq "TGT" or $_[0] eq "TGC");
	return "A" if($_[0] eq "GCT" or $_[0] eq "GCC" or $_[0] eq "GCA" or $_[0] eq "GCG");
	return "G" if($_[0] eq "GGT" or $_[0] eq "GGC" or $_[0] eq "GGA" or $_[0] eq "GGG");
	return "P" if($_[0] eq "CCT" or $_[0] eq "CCC" or $_[0] eq "CCA" or $_[0] eq "CCG");
	return "T" if($_[0] eq "ACT" or $_[0] eq "ACC" or $_[0] eq "ACA" or $_[0] eq "ACG");
	return "S" if($_[0] eq "TCT" or $_[0] eq "TCC" or $_[0] eq "TCA" or $_[0] eq "TCG" or $_[0] eq "AGT" or $_[0] eq "AGC");
	return "Y" if($_[0] eq "TAT" or $_[0] eq "TAC");
	return "W" if($_[0] eq "TGG");
	return "Q" if($_[0] eq "CAA" or $_[0] eq "CAG");
	return "N" if($_[0] eq "AAT" or $_[0] eq "AAC");
	return "H" if($_[0] eq "CAT" or $_[0] eq "CAC");
	return "E" if($_[0] eq "GAA" or $_[0] eq "GAG");
	return "D" if($_[0] eq "GAT" or $_[0] eq "GAC");
	return "K" if($_[0] eq "AAA" or $_[0] eq "AAG");
	return "R" if($_[0] eq "CGT" or $_[0] eq "CGC" or $_[0] eq "CGA" or $_[0] eq "CGG" or $_[0] eq "AGA" or $_[0] eq "AGG");
	return "*" if($_[0] eq "TAA" or $_[0] eq "TAG" or $_[0] eq "TGA");
}

sub n_of_codon
{
    return 3 if($_[0] eq "ATT" or $_[0] eq "ATC" or $_[0] eq "ATA");
    return 6 if($_[0] eq "CTT" or $_[0] eq "CTC" or $_[0] eq "CTA" or $_[0] eq "CTG" or $_[0] eq "TTA" or $_[0] eq "TTG");
    return 4 if($_[0] eq "GTT" or $_[0] eq "GTC" or $_[0] eq "GTA" or $_[0] eq "GTG");
    return 2 if($_[0] eq "TTT" or $_[0] eq "TTC");
    return 1 if($_[0] eq "ATG");
    return 2 if($_[0] eq "TGT" or $_[0] eq "TGC");
    return 4 if($_[0] eq "GCT" or $_[0] eq "GCC" or $_[0] eq "GCA" or $_[0] eq "GCG");
    return 4 if($_[0] eq "GGT" or $_[0] eq "GGC" or $_[0] eq "GGA" or $_[0] eq "GGG");
    return 4 if($_[0] eq "CCT" or $_[0] eq "CCC" or $_[0] eq "CCA" or $_[0] eq "CCG");
    return 4 if($_[0] eq "ACT" or $_[0] eq "ACC" or $_[0] eq "ACA" or $_[0] eq "ACG");
    return 6 if($_[0] eq "TCT" or $_[0] eq "TCC" or $_[0] eq "TCA" or $_[0] eq "TCG" or $_[0] eq "AGT" or $_[0] eq "AGC");
    return 2 if($_[0] eq "TAT" or $_[0] eq "TAC");
    return 1 if($_[0] eq "TGG");
    return 2 if($_[0] eq "CAA" or $_[0] eq "CAG");
    return 2 if($_[0] eq "AAT" or $_[0] eq "AAC");
    return 2 if($_[0] eq "CAT" or $_[0] eq "CAC");
    return 2 if($_[0] eq "GAA" or $_[0] eq "GAG");
    return 2 if($_[0] eq "GAT" or $_[0] eq "GAC");
    return 2 if($_[0] eq "AAA" or $_[0] eq "AAG");
    return 6 if($_[0] eq "CGT" or $_[0] eq "CGC" or $_[0] eq "CGA" or $_[0] eq "CGG" or $_[0] eq "AGA" or $_[0] eq "AGG");
    return 3 if($_[0] eq "TAA" or $_[0] eq "TAG" or $_[0] eq "TGA");
}


#input sequence must be of the same length, no exception is concerned, and the length of which is divisible by 3.
#I need to count each codon(hash) and each pep in order to
#calculate RSCU
#I need to count each codon(hash) and each pep in order to
#calculate RSCU
sub calc_GC3_RSCU
{
	my ($count_codon, $count_gc3, $GC3) = (0, 0, 0);

	undef %RSCU;
	$GC3="E";

	$input1=shift;
	
	@seq1_nuc=split //, $$input1;

	$len_1_3=int(@seq1_nuc/3);
	
	for($i_1_3=0; $i_1_3<$len_1_3; $i_1_3++)
	{
		$i_pos1=$i_1_3*3;
		$i_pos2=$i_pos1+1;
		$i_pos3=$i_pos2+1;

		$codon1=$seq1_nuc[$i_pos1].$seq1_nuc[$i_pos2].$seq1_nuc[$i_pos3];

		$codon1 =~ s/(.*)/\U$1/g;

		next if($codon1=~/N/);
		
		$pep1=&nuc2pep($codon1);

		$Counter{$codon1}++;
		$Counter{$pep1}++;
		$Counter_Total{$codon1}++;
		$Counter_Total{$pep1}++;

		$count_codon++;
		$count_gc3++ if ($seq1_nuc[$i_pos3] eq "G" or $seq1_nuc[$i_pos3] eq "C");

		$Counter_GC3_Total++ if ($seq1_nuc[$i_pos3] eq "G" or $seq1_nuc[$i_pos3] eq "C");
		$Counter_Codon_Total++;

	}
	$GC3=$count_gc3/$count_codon if ($count_codon !=0);
	for($i_name=0; $i_name<64; $i_name++)
	{
		if(defined($Counter{&nuc2pep($codon_name[$i_name])}))
		{
			if(defined($Counter{$codon_name[$i_name]}))
			{
				$RSCU{$codon_name[$i_name]}=$Counter{$codon_name[$i_name]}/
					$Counter{&nuc2pep($codon_name[$i_name])}*&n_of_codon($codon_name[$i_name]);
			}
			else
			{
				$RSCU{$codon_name[$i_name]}=0;
			}
		}
		else
		{
			$RSCU{$codon_name[$i_name]}=1;
		}
	}
	return $GC3;

}

#func: print gene name GC3 64 codon's RSCU
sub calc_all
{
	$gc3=&calc_GC3_RSCU(\$sequence);
	print "$mem_gene_name\t$gc3\t";
#	for($i_name=0; $i_name<64; $i_name++)
#	{
#		print $RSCU{$codon_name[$i_name]}."\t";
#	}
	for $aa_tmp(sort keys %aminoacid)
	{
		for $codon_tmp (@{$aminoacid{$aa_tmp}})
		{
			printf ("%.2f\t",  $RSCU{$codon_name[$number_codon{$codon_tmp}]});
		}
	}
	print "\n";
}

#inherited from sub calc_all
sub calc_using_xx_Total
{
	$GC3_Total = $Counter_GC3_Total/$Counter_Codon_Total;
        for($i_name=0; $i_name<64; $i_name++)
        {
                if(defined($Counter_Total{&nuc2pep($codon_name[$i_name])}))
                {
                        if(defined($Counter_Total{$codon_name[$i_name]}))
                        {
                                $RSCU_Total{$codon_name[$i_name]}=$Counter_Total{$codon_name[$i_name]}/
                                $Counter_Total{&nuc2pep($codon_name[$i_name])}*&n_of_codon($codon_name[$i_name]);
                        }
                        else
                        {
                                $RSCU_Total{$codon_name[$i_name]}=0;
                        }
                }
                else
                {
                        $RSCU_Total{$codon_name[$i_name]}=1;
                }
        }

        print "Total\t$GC3_Total\t";

        for $aa_tmp(sort keys %aminoacid)
        {
                for $codon_tmp (@{$aminoacid{$aa_tmp}})
                {
                        printf ("%.2f\t",  $RSCU_Total{$codon_name[$number_codon{$codon_tmp}]});
                }
        }
        print "\n";
}


#sub is over, main program began
#print "                   Indica_Nivara\tJaponica_Rufipogon\n";
#print "Gene        \tN_S\tS\tSNP\tN_S\tS\tSNP\n";

$_=<>;
@e=split;
$mem_gene_name=substr($e[0],1);
$sequence="";


###################################################################
###################################################################
#print header
print "name\tGC3\t";
for $aa_tmp(sort keys %aminoacid)
{
        for $codon_tmp (@{$aminoacid{$aa_tmp}})
        {
                print $aa_tmp."\t";
        }
}
print "\n";
#print sec-header
print "name\tGC3\t";
for $aa_tmp(sort keys %aminoacid)
{
       for $codon_tmp (@{$aminoacid{$aa_tmp}})
        {
		print $codon_tmp."\t";
	}
}
print "\n";
##################################################################
##################################################################



#read gene sequence, 50 species, genes with gaps are also read into memory
while(<>)
{
	chomp;
	if($_=~/>/)
	{
		$next_gene_name=substr($_, 1);
		&calc_all;#will use $mem_gene_name and sequence as hidden param
		$mem_gene_name=$next_gene_name;
		
		$sequence="";
	}
	else
	{
		s/(.*)/$1\U/g;
		$sequence.=$_;
	}
}

&calc_all;#will use $mem_gene_name and sequence as hidden param
&calc_using_xx_Total;


