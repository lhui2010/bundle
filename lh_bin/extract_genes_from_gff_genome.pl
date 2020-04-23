#!/usr/bin/perl

#version 1.0 extract genes from gff's loci(generated by ~/.../TOTAL/toloci.pl) and genome's fasta
#func by reading loci data

#fixed a bug in rc subroutine


if(@ARGV < 2)
{
	print "usage:$0 ChrXX'sgff'sLoci ChrXX'sgenome >gene\n" if(@ARGV == 0);
	exit;
}

open Floci, $ARGV[0] or die$!;
open Fgenome, $ARGV[1] or die $!;


$genome = ">";

$info_fa=<Fgenome>;
while(<Fgenome>)
{
	chomp;
	$genome.=$_;
}
close Fgenome;

while(<Floci>)
{
	my (@e, $chr, $strand, $genename, $cdsseq);
        chomp;
	@e=split;
	$chr = shift @e;
	$strand = shift @e;
	$genename = shift @e;

	for $loci(@e)
	{
		
		my ($start, $end)=split /-/, $loci;
		my $exonlen = $end - $start +1;
		my $exonseq = substr($genome, $start, $exonlen);
		
		if($strand eq "-")
		{
			&rc (\$exonseq);
		}

		$cdsseq .= $exonseq;
	}
	
	print ">$genename\n";
	&formated_print (\$cdsseq);
}
close Floci;







#----SUBROUTINE----


sub formated_print
{
	        $length = 70;
		        my $seq_p =shift;
			        if (ref($seq_p) eq 'SCALAR') {
					                for ( my $pos = 0 ; $pos < length($$seq_p) ; $pos += $length ) {
								                        print substr($$seq_p, $pos, $length), "\n";
											                }
													        }
}

#sub Complement_Reverse
sub rc
{
	        my $seq_p=shift;
		        if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
				                $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
						                $$seq_p=reverse($$seq_p);
								        }
}


