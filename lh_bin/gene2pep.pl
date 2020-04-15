#!/usr/bin/perl -w

#genes have been selected, it's calculation time

print "usage: .pl gene.fa >pep.fa
output: pep.fa
"if(@ARGV==0);

while(<>)
{
	chomp;
	if(/>/)
	{
		if(defined($gene_name))
		{
			print '>',$gene_name,"\n";
			&formated_print(&nuc2pep_seq($seq));
			$seq="";
		}
		$gene_name=substr($_,1);
	}
	else
	{
#		$_="aaaaaaaaTTTaaaaaaa";
#		print $_,"\n";
		s/$_/\U$_/gi;#print $_,"\n";exit;
		$seq.=$_;
	}
}

print '>',$gene_name,"\n";
&formated_print(&nuc2pep_seq($seq));

#------------sub------------

sub formated_print
{
    $len_input=length($_[0]);
    for ($no=0; $no<$len_input; $no+=70)
    {
        print substr($_[0], $no, 70),"\n";
    }
}


sub nuc2pep_seq
{
	my ($pep_seq, $len);
	$len=length($_[0]);
	for($i=0; $i<$len; $i+=3)
	{
		$pep_seq.=&nuc2pep(substr($_[0], $i, 3));
	}
	return $pep_seq;
}

sub count_codon_per_aa
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
	return "-" if($_[0]=~/-/);
}

	
