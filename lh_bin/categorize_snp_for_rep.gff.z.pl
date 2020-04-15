#!/usr/bin/perl

#modified by liuhui to adjust .gz input
#Sun Apr 17 16:15:03 CST 2011

=head1 Name

rice_reseq_snp_distri.pl -- detect the snp distribution in the genome: intergenic region, genic region, 5'-UTR, 3'-URT, intron,CDS, (synonymous and nonsynonymous)

=head1 Description 

read the whole population snp data, gff data, genome reference, and generate the files which contain different distribution

=head1 Version
	
	Author: Lv Jun, lvu.jun@gmail.com	(mimic Fan Wei)
	Version: 1.0, Date: 2009-10-03
	Note: for rice

=head1 Usage

	--verbose output running progress information to screen
	--help	  output help information to screen

=head1 example

	perl rice_reseq_snp_distri.pl /home/lvjun/data/new/IRGSP_ref/analysis/no_LC_filter/whole/chromosome01 /home/lvjun/data/GFF/gff_RAP2/rep.gff /home/lvjun/data/30Rice_Resequencing/reference/IRGSP_chr/Chr01 -o ./output

=cut

use strict; 
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;  ## function "mkpath" and "rmtree" deal with directory

## get options form command line into variables and set default values
my ($Verbose,$Help,$Outdir);
GetOptions(
		"outdir:s"=>\$Outdir,
		"verbose"=>\$Verbose,
		"help"=>\$Help

);
$Outdir ||= "./";
die `pod2text $0` if (@ARGV == 0 || $Help);

## to read in the filenames from the command line
my $base_file = shift; 
my $gff_file = shift;
my $genome_ref = shift;


$Outdir =~ s/\/$//;  # if the last sign is /, then delete it
my $base_file_basename = basename($base_file);
my $genic_output_info_file = "$Outdir/$base_file_basename.genic.info";
#my $genic_output_seq_file = "$Outdir/$base_file_basename_genic.seq";
#my $genic_output_axt_file = "$Outdir/$base_file_basename_genic.axt";
my $UTR5_output_info_file = "$Outdir/$base_file_basename.UTR5.info";
#my $UTR5_output_seq_file = "$Outdir/$base_file_basename_UTR5.seq";
#my $UTR5_output_axt_file = "$Outdir/$base_file_basename_UTR5.axt";
my $UTR3_output_info_file = "$Outdir/$base_file_basename.UTR3.info";
my $exon_output_info_file = "$Outdir/$base_file_basename.exon.info";
#my $UTR3_output_seq_file = "$Outdir/$base_file_basename_UTR3.seq";
#my $UTR3_output_axt_file = "$Outdir/$base_file_basename_UTR3.axt";
my $mRNA_output_info_file = "$Outdir/$base_file_basename.mRNA.info";
my $cds_output_info_file = "$Outdir/$base_file_basename.cds.info";
my $cds_output_seq_file = "$Outdir/$base_file_basename.cds.seq";
my $cds_output_axt_file = "$Outdir/$base_file_basename.cds.axt";

#debug:
#print $_, "\n" for ($base_file_basename, $genic_output_info_file, $UTR5_output_info_file, $UTR3_output_info_file, $exon_output_info_file, $mRNA_output_info_file, $cds_output_info_file, $cds_output_seq_file, $cds_output_axt_file);
#exit;

my %SNP; ## to store SNP data
my %gene; ## store gene data
my %CDS; ## store cds data
my %UTR5; ## store 5-UTR data
my %UTR3; ## store 3-UTR data
my %mRNA;
my %exon;

my %CODE = (
	'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
	'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
	'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
	'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
	'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
	'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
	'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
	'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
	'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
	'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
	'ATG' => 'M',                                                                         # Methionine
	'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
	'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
	'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
	'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
	'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
	'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
	'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
	'TGG' => 'W',                                                                         # Tryptophan
	'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
	'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
);

my %Abbrev = (
		'A' => [ 'A' ],
		'C' => [ 'C' ],
		'G' => [ 'G' ],
		'T' => [ 'T' ],
		'M' => [ 'A', 'C' ],
		'R' => [ 'A', 'G' ],
		'W' => [ 'A', 'T' ],
		'S' => [ 'C', 'G' ],
		'Y' => [ 'C', 'T' ],
		'K' => [ 'G', 'T' ],
		'V' => [ 'A', 'C', 'G' ],
		'H' => [ 'A', 'C', 'T' ],
		'D' => [ 'A', 'G', 'T' ],
		'B' => [ 'C', 'G', 'T' ],
		'X' => [ 'A', 'C', 'G', 'T' ],  
		'N' => [ 'A', 'C', 'G', 'T' ]
);

## chromosome01    51      T       4543    3.453115        T       A       49      1       3       16.68   1587    3285 
my %capui = ("AC","M", "CA","M", "GT","K", "TG","K", "CT","Y", "TC","Y", "AG","R", "GA","R", "AT","W", "TA","W", "CG","S", "GC","S");

open SNP, "gzip -dc $base_file|" || die "can't open the $base_file\n";
while(<SNP>)
{
	chomp;
	my @temp = split;
	my $chr = $temp[0];
	my $pos = $temp[1];
	my $ref_type = $temp[2];
	my $rice_type = $temp[3];
	$SNP{$chr}{$pos} = [$ref_type,$rice_type];
		
}
close(SNP);

warn "read base_file done" if($Verbose);

#chr1    mRNA    gene    287362  291373  .       +       .       locus_id "Os01g0105400";Description "Similar to Kinesin heavy chain.";category "II : Similar
#chr1    mRNA    rep_CDS 24883   25391   .       +       1       seq_id "AK121362";locus_id "Os01g0100600"
#chr1    mRNA    rep_5'-UTR      26143   26220   .       +       .       seq_id "AK121523";locus_id "Os01g0100700"
#chr1    mRNA    rep_3'-UTR      27420   27644   .       +       .       seq_id "AK121523";locus_id "Os01g0100700"
#chr1    COMBINER_EST    rep_EST_predicted_CDS   351390  351459  .       +       2       seq_id "Os01g0106700";locus_id "Os01g0106700"
#chr1    mRNA    rep_UTR 8631904 8632355 .       -       .       seq_id "AK101657";locus_id "Os01g0258100";type "Internal"

open (GFF,$gff_file) || die "can't open the $gff_file\n";
while(<GFF>)
{
	chomp;
	my @temp = split(/\t/);
#	my $temp = join "\t",@temp[0..$#temp];
	my $chr_temp = $temp[0];
	my $num = $1 if($chr_temp =~/(\d+)/);
	my $chr;
	$chr = "chromosome$num";
	my $type = $temp[2];
	my $start = $temp[3];
	my $end = $temp[4];
	my $strand = $temp[6];
	my $ID;
#	print "$type\t$start\t$end\n";
#	print "$temp\n$type\n";
		if($type =~ /gene/) { $ID = $1 if($temp[8] =~/^ID=(\S+);Name/); push @{$gene{$chr}{$ID}},[$start,$end,$strand];}
		elsif($type =~/CDS/) { $ID = $1 if ($temp[8] =~/Parent=(\S+)/);push @{$CDS{$chr}{$ID}},[$start,$end,$strand];}
		elsif($type =~/exon/) { $ID = $1 if ($temp[8] =~/Parent=(\S+)/);push @{$exon{$chr}{$ID}},[$start,$end,$strand];}
		elsif($type =~/five_prime_UTR/) { $ID = $1 if ($temp[8] =~/Parent=(\S+)/);push @{$UTR5{$chr}{$ID}},[$start,$end,$strand];}
		elsif($type =~/three_prime_UTR/) { $ID = $1 if ($temp[8] =~/Parent=(\S+)/);push @{$UTR3{$chr}{$ID}},[$start,$end,$strand];}
		elsif($type =~/mRNA/) {$ID = $1 if($temp[8] =~/ID=(\S+);Name/);push @{$mRNA{$chr}{$ID}},[$start,$end,$strand];}

}
close(GFF);
warn "read GFF done" if($Verbose);
#foreach (keys %{$pre_CDS{chromosome01}}) {print "$_\n";}
#exit;

open INFO_GENE, ">$genic_output_info_file" or die "fail output to $genic_output_info_file";
open INFO_UTR5, ">$UTR5_output_info_file" or die "fail output to $UTR5_output_info_file";
open INFO_UTR3, ">$UTR3_output_info_file" or die "fail output to $UTR3_output_info_file";
open INFO_EXON, ">$exon_output_info_file" or die "fail output to $exon_output_info_file";
open INFO_mRNA, ">$mRNA_output_info_file" or die "fail output to $mRNA_output_info_file";
open INFO_CDS, ">$cds_output_info_file" or die "fail output to $cds_output_info_file";
open SEQ_CDS, ">$cds_output_seq_file" or die "fail output to $cds_output_seq_file";
open AXT_CDS, ">$cds_output_axt_file" or die "fail output to $cds_output_axt_file";

print INFO_CDS "chr_id\tposition\tref_base<->reseq_base\tsnp_status\tstrand\tAnnotation_type\tmRNA_id\tGene_id\tFeature_type\tcodon_phase\tcodon_mutate\taa_mutate\tsynonymous\tnonsynonymous\n";
print INFO_GENE "chr_id\tposition\tref_base<->reseq_base\tsnp_status\tstrand\tAnnotation_type\tGene_id\tFeature_type\n";
print INFO_UTR5 "chr_id\tposition\tref_base<->reseq_base\tsnp_status\tstrand\tAnnotation_type\tmRNA_id\tGene_id\tFeature_type\n";
print INFO_UTR3 "chr_id\tposition\tref_base<->reseq_base\tsnp_status\tstrand\tAnnotation_type\tmRNA_id\tGene_id\tFeature_type\n";
print INFO_EXON "chr_id\tposition\tref_base<->reseq_base\tsnp_status\tstrand\tAnnotation_type\tmRNA_id\tGene_id\tFeature_type\n";
print INFO_mRNA "chr_id\tposition\tref_base<->reseq_base\tsnp_status\tstrand\tAnnotation_type\tmRNA_id\tGene_id\tFeature_type\n";

##read in the reference genome sequence(eg. IRGSP)
open REF, $genome_ref || die "can't open $genome_ref \n";
$/ = ">"; <REF>; $/ = "\n";
while(<REF>)
{
	my $chr = $1 if(/^(\S+)/);
	$/ = ">";
	my $seq = <REF>;
	chomp $seq;
	$seq =~s/\s//g;
	$seq = uc($seq); ## all upper case
	$/ = "\n";

	warn "read rice genome $chr done"if ($Verbose);  ## to deal with the chromosome which has just been read in
	next if (!exists $SNP{$chr});
	
	my $chr_snp_p = $SNP{$chr};
	my $chr_gene_p = $gene{$chr};
	my $chr_cds_p = $CDS{$chr};
	my $chr_utr5 = $UTR5{$chr};
	my $chr_utr3 = $UTR3{$chr};
	my $chr_exon = $exon{$chr};
	my $chr_mRNA_p = $mRNA{$chr};

	my ($output_info_gene,$output_info_utr5,$output_info_utr3,$output_info_exon,$output_info_cds,$output_axt_cds,$output_seq_cds,$output_info_mRNA);

	foreach my $gene_id (keys %$chr_gene_p)			### to deal with every single gene on this chromosome
		{
			my $gene_p = $chr_gene_p->{$gene_id};
			my $strand = $gene_p->[0][2];			### actually for every gene, there's only one record, unlike CDS
		     for( my $j = 0 ; $j <@$gene_p; $j++)		# in case that the same gene name can have several interval, like alternative splicing
		      	{
			for (my $i = $gene_p->[$j][0];$i<=$gene_p->[$j][1]; $i++)			# loop for every base
				{	
					my $ref_base;
					if(exists $chr_snp_p->{$i}) 
					{ $ref_base = substr($seq,$i-1,1);
					warn "ref_base error on $chr:$i, please check: $ref_base $chr_snp_p->{$i}[0]" if($ref_base ne $chr_snp_p->{$i}[0]);

					## change to the + strand and then calculate
					my $rice_base = $chr_snp_p->{$i}[1];	# the [0] is the reference type, [1] is the reseq type
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_rice_base = ($strand eq "+") ? $rice_base : $rice_base_minus;
					my $this_ref_base = ($strand eq "+") ? $ref_base : $ref_base_minus;	### $ref_base belongs to + strand
					my $snp_status = hom_het($this_ref_base,$this_rice_base);
					$output_info_gene .= "$chr\t$i\t$ref_base<->$rice_base\t$snp_status\t$strand\tgene\t$gene_id\tgene\n";
					}
				}
			}
		}	
		print INFO_GENE $output_info_gene;

	foreach my $gene_id (keys %$chr_utr5)			### to deal with every single utr5 on this chromosome
		{
			my $utr5_p = $chr_utr5->{$gene_id};
			my $strand = $utr5_p->[0][2];			
		  for (my $j = 0; $j < @$utr5_p; $j++)			# for one gene, there can be several utr5 segments
		       {
			for (my $i = $utr5_p->[$j][0];$i<=$utr5_p->[$j][1]; $i++)			# loop for every base
				{	
					my $ref_base;
					if(exists $chr_snp_p->{$i}) 
					{ $ref_base = substr($seq,$i-1,1);
					warn "ref_base error on $chr:$i, please check: $ref_base $chr_snp_p->{$i}[0]" if($ref_base ne $chr_snp_p->{$i}[0]);

					## change to the + strand and then calculate
					my $rice_base = $chr_snp_p->{$i}[1];	# the [0] is the reference type, [1] is the reseq type
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_rice_base = ($strand eq "+") ? $rice_base : $rice_base_minus;
					my $this_ref_base = ($strand eq "+") ? $ref_base : $ref_base_minus;
					my $snp_status = hom_het($this_ref_base,$this_rice_base);
					$output_info_utr5 .= "$chr\t$i\t$ref_base<->$rice_base\t$snp_status\t$strand\tgene\t$gene_id\tutr5\n";
					}
				}
			}
		}

		print INFO_UTR5 $output_info_utr5;

	foreach my $gene_id (keys %$chr_utr3)			### to deal with every single utr3 on this chromosome
		{
			my $utr3_p = $chr_utr3->{$gene_id};
			my $strand = $utr3_p->[0][2];			
		for (my $j = 0; $j < @$utr3_p;$j++)			### for one gene, there can be several utr5 segments
		       {
			for (my $i = $utr3_p->[$j][0];$i<=$utr3_p->[$j][1]; $i++)			# loop for every base
				{	
					my $ref_base;
					if(exists $chr_snp_p->{$i}) 
					{ $ref_base = substr($seq,$i-1,1);
					warn "ref_base error on $chr:$i, please check: $ref_base $chr_snp_p->{$i}[0]" if($ref_base ne $chr_snp_p->{$i}[0]);

					## change to the + strand and then calculate
					my $rice_base = $chr_snp_p->{$i}[1];	# the [0] is the reference type, [1] is the reseq type
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_rice_base = ($strand eq "+") ? $rice_base : $rice_base_minus;
					my $this_ref_base = ($strand eq "+") ? $ref_base : $ref_base_minus;
					my $snp_status = hom_het($this_ref_base,$this_rice_base);
					$output_info_utr3 .= "$chr\t$i\t$ref_base<->$rice_base\t$snp_status\t$strand\tgene\t$gene_id\tutr3\n";
					}
				}
			}
		}
			
	print INFO_UTR3 $output_info_utr3;

	foreach my $gene_id (keys %$chr_exon)			### to deal with every single exon on this chromosome
		{
			my $exon_p = $chr_exon->{$gene_id};
			my $strand = $exon_p->[0][2];			
		for (my $j = 0; $j < @$exon_p;$j++)			### for one gene, there can be several utr5 segments
		       {
			for (my $i = $exon_p->[$j][0];$i<=$exon_p->[$j][1]; $i++)			# loop for every base
				{	
					my $ref_base;
					if(exists $chr_snp_p->{$i}) 
					{ $ref_base = substr($seq,$i-1,1);
					warn "ref_base error on $chr:$i, please check: $ref_base $chr_snp_p->{$i}[0]" if($ref_base ne $chr_snp_p->{$i}[0]);

					## change to the + strand and then calculate
					my $rice_base = $chr_snp_p->{$i}[1];	# the [0] is the reference type, [1] is the reseq type
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_rice_base = ($strand eq "+") ? $rice_base : $rice_base_minus;
					my $this_ref_base = ($strand eq "+") ? $ref_base : $ref_base_minus;
					my $snp_status = hom_het($this_ref_base,$this_rice_base);
					$output_info_exon .= "$chr\t$i\t$ref_base<->$rice_base\t$snp_status\t$strand\tgene\t$gene_id\texon\n";
					}
				}
			}
		}
			
	print INFO_EXON $output_info_exon;


	foreach my $gene_id (keys %$chr_mRNA_p)			### to deal with every single mRNA on this chromosome
		{
			my $mRNA_p = $chr_mRNA_p->{$gene_id};
			my $strand = $mRNA_p->[0][2];			### actually for every gene, there's only one record, unlike CDS
		  for( my $j = 0; $j < @$mRNA_p; $j++)
		       {	
			for (my $i = $mRNA_p->[$j][0];$i<=$mRNA_p->[$j][1]; $i++)			# loop for every base
				{	
					my $ref_base;	
					if(exists $chr_snp_p->{$i}) 
					{$ref_base = substr($seq,$i-1,1);
					warn "ref_base error on $chr:$i, please check: $ref_base $chr_snp_p->{$i}[0]" if($ref_base ne $chr_snp_p->{$i}[0]);

					## change to the + strand and then calculate
					my $rice_base = $chr_snp_p->{$i}[1];	# the [0] is the reference type, [1] is the reseq type
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_rice_base = ($strand eq "+") ? $rice_base : $rice_base_minus;
					my $this_ref_base = ($strand eq "+") ? $ref_base : $ref_base_minus;
					my $snp_status = hom_het($this_ref_base,$this_rice_base);
					$output_info_mRNA .= "$chr\t$i\t$ref_base<->$rice_base\t$snp_status\t$strand\tgene\t$gene_id\tmRNA\n";
					}
				}
			}
		}
	print INFO_mRNA $output_info_mRNA;


	foreach my $gene_id (keys %$chr_cds_p) 			## to calculate the CDS, syn and non-syn SNP
	{
		my $gene_p = $chr_cds_p->{$gene_id};
		my $strand = $gene_p->[0][2];
		my $gene_cds_str;
		my $gene_cds_len = 0;
		my $mRNA_pos = 1;
		my $is_cds_mutated = 0;
		my $ref_cds_str;
		my $axt_cds_str; ##if meet het-two SNP, only keep one by random
		
		##.....cds.............................
	    if ($strand eq "+")
	      {
		for (my $i=0; $i<@$gene_p; $i++) 
		{ ##loop for each exon
			my $cds_p = $gene_p->[$i];
			$gene_cds_str .= substr($seq,$cds_p->[0]-1,$cds_p->[1]-$cds_p->[0]+1);
		}
	      }
	    elsif ($strand eq "-")
	      {
		for (my $i=@$gene_p-1; $i >= 0; $i--)
		{ ##loop for each exon 
			my $cds_p = $gene_p->[$i];
			$gene_cds_str .= substr($seq,$cds_p->[0]-1,$cds_p->[1]-$cds_p->[0]+1);
		}
		
	      }
		$gene_cds_len = length($gene_cds_str);
		Complement_Reverse(\$gene_cds_str) if ($strand eq "-");
		$ref_cds_str = $gene_cds_str;
		$axt_cds_str = $gene_cds_str;
	    if ($strand eq "+")
	      {
		for (my $i=0; $i<@$gene_p; $i++) 
		{ ##loop for each exon
			my $cds_p = $gene_p->[$i];
			for (my $j=$cds_p->[0]; $j<=$cds_p->[1]; $j++) 
			{ ##loop for each base
				if (exists $chr_snp_p->{$j})
				{
					$is_cds_mutated = 1;
					my $ref_base = substr($seq,$j-1,1);
					warn "ref_base error on $chr:$j, please check: $ref_base $chr_snp_p->{$j}[0]" if($ref_base ne $chr_snp_p->{$j}[0]);
					
					##................
					my $rice_base = $chr_snp_p->{$j}[1];
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_mRNA_pos =  $mRNA_pos ;
					my $this_rice_base =  $rice_base ;
					my $this_ref_base =  $ref_base ;

					substr($gene_cds_str,$this_mRNA_pos-1,1) = $this_rice_base;
#print"this_rice_base: $this_rice_base";exit;
					substr($axt_cds_str,$this_mRNA_pos-1,1) = different_base($this_ref_base,$this_rice_base);
					my $snp_status = hom_het($ref_base,$rice_base);
					
					##chr10   104221044    A<->R   Het-one + Gene    NM_024789        CDS    729:2       GTA<->GTG;      V<->V;  0.5     0
					$output_info_cds .= "$chr\t$j\t$ref_base<->$rice_base\t$snp_status\t$strand\tGene\t";
					
					my ($codon_num,$phase_num) = codon_phase($this_mRNA_pos);
					my $codon_phase_str = "$this_mRNA_pos:$phase_num";
					my $ref_codon = substr($ref_cds_str,$this_mRNA_pos-$phase_num-1,3);
					my $rice_codon = substr($gene_cds_str,$this_mRNA_pos-$phase_num-1,3);
					
					warn "$rice_codon has N, please check\n" if($rice_codon =~ /N/);
					
					##get the yanhuang two codons
					my ($rice_codon1,$rice_codon2) = ("","");
					($rice_codon1,$rice_codon2)= convert_codon($rice_codon,$ref_codon,$phase_num) if($snp_status =~ /Het/);
					$rice_codon1 = $rice_codon2 = $rice_codon if($snp_status =~ /Hom/);
					
					##caculate the synonymous and nonsynonymous SNP number
					my ($synonymous,$nonsynonymous) = (0,0);
					my ($codon_mutate_str,$aa_mutate_str);
					if ($ref_codon ne $rice_codon1) 
					{
						if ($CODE{$ref_codon} eq $CODE{$rice_codon1}) {
							$synonymous += 0.5;
						}else{
							$nonsynonymous += 0.5;
						}
						$codon_mutate_str .= "$ref_codon<->$rice_codon1;";
						$aa_mutate_str .= "$CODE{$ref_codon}<->$CODE{$rice_codon1};";
					}
					if ($ref_codon ne $rice_codon2) 
					{
						if ($CODE{$ref_codon} eq $CODE{$rice_codon2}) {
							$synonymous += 0.5;
						}else{
							$nonsynonymous += 0.5;
						}
						$codon_mutate_str .= "$ref_codon<->$rice_codon2;" if($rice_codon2 ne $rice_codon1);
						$aa_mutate_str .= "$CODE{$ref_codon}<->$CODE{$rice_codon2};" if($rice_codon2 ne $rice_codon1);
					}

					$output_info_cds .= "$gene_id\tCDS\t$codon_phase_str\t$codon_mutate_str\t$aa_mutate_str\t$synonymous\t$nonsynonymous\n";
				}

				$mRNA_pos++; ##postion on strand + 

			}			
		}
	   }
	
	else 	
		 {
		for (my $i= @$gene_p -1; $i >= 0; $i--) 
		{ ##loop for each exon
			my $cds_p = $gene_p->[$i];
			for (my $j=$cds_p->[0]; $j<=$cds_p->[1]; $j++) 
			{ ##loop for each base
				if (exists $chr_snp_p->{$j})
				{
					$is_cds_mutated = 1;
					my $ref_base = substr($seq,$j-1,1);
					warn "ref_base error on $chr:$j, please check: $ref_base $chr_snp_p->{$j}[0]" if($ref_base ne $chr_snp_p->{$j}[0]);
					
					##................
					my $rice_base = $chr_snp_p->{$j}[1];
					my $rice_base_minus = $rice_base;
					$rice_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $ref_base_minus = $ref_base;
					$ref_base_minus =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					my $this_mRNA_pos = $gene_cds_len - $mRNA_pos + 1 ;
					my $this_rice_base = $rice_base_minus;
					my $this_ref_base = $ref_base_minus;

					substr($gene_cds_str,$this_mRNA_pos-1,1) = $this_rice_base;
#print"this_rice_base: $this_rice_base";exit;
					substr($axt_cds_str,$this_mRNA_pos-1,1) = different_base($this_ref_base,$this_rice_base);
					my $snp_status = hom_het($ref_base,$rice_base);
					
					##chr10   104221044    A<->R   Het-one + Gene    NM_024789        CDS    729:2       GTA<->GTG;      V<->V;  0.5     0
					$output_info_cds .= "$chr\t$j\t$ref_base<->$rice_base\t$snp_status\t$strand\tGene\t";
					
					my ($codon_num,$phase_num) = codon_phase($this_mRNA_pos);
					my $codon_phase_str = "$this_mRNA_pos:$phase_num";
					my $ref_codon = substr($ref_cds_str,$this_mRNA_pos-$phase_num-1,3);
					my $rice_codon = substr($gene_cds_str,$this_mRNA_pos-$phase_num-1,3);
					
					warn "$rice_codon has N, please check\n" if($rice_codon =~ /N/);
					
					##get the yanhuang two codons
					my ($rice_codon1,$rice_codon2);
					($rice_codon1,$rice_codon2)= convert_codon($rice_codon,$ref_codon,$phase_num) if($snp_status =~ /Het/);
					$rice_codon1 = $rice_codon2 = $rice_codon if($snp_status =~ /Hom/);
					
					##caculate the synonymous and nonsynonymous SNP number
					my ($synonymous,$nonsynonymous) = (0,0);
					my ($codon_mutate_str,$aa_mutate_str);
					if ($ref_codon ne $rice_codon1) 
					{
						if ($CODE{$ref_codon} eq $CODE{$rice_codon1}) {
							$synonymous += 0.5;
						}else{
							$nonsynonymous += 0.5;
						}
						$codon_mutate_str .= "$ref_codon<->$rice_codon1;";
						$aa_mutate_str .= "$CODE{$ref_codon}<->$CODE{$rice_codon1};";
					}
					if ($ref_codon ne $rice_codon2) 
					{
						if ($CODE{$ref_codon} eq $CODE{$rice_codon2}) {
							$synonymous += 0.5;
						}else{
							$nonsynonymous += 0.5;
						}
						$codon_mutate_str .= "$ref_codon<->$rice_codon2;" if($rice_codon2 ne $rice_codon1);
						$aa_mutate_str .= "$CODE{$ref_codon}<->$CODE{$rice_codon2};" if($rice_codon2 ne $rice_codon1);
					}

					$output_info_cds .= "$gene_id\tCDS\t$codon_phase_str\t$codon_mutate_str\t$aa_mutate_str\t$synonymous\t$nonsynonymous\n";
				}

				$mRNA_pos++; ##postion on strand + 

			}			
		}
	   }
	
##.........is_matated...........
		$output_axt_cds .="$gene_id\_ref&$gene_id\_rice\n".$ref_cds_str."\n".$axt_cds_str."\n\n" ;

		Display_seq(\$ref_cds_str);
#		$output_refcds .= ">$gene_id\n".$ref_cds_str;
		
		Display_seq(\$gene_cds_str);
		$output_seq_cds .= ">$gene_id\n".$gene_cds_str;
		##warn "SNP caculation for $gene_id done" if($Verbose);
	}
	
	warn "SNP caculation for $chr done" if($Verbose);
	print INFO_CDS $output_info_cds;
#	print REFCDS $output_refcds;
	print SEQ_CDS $output_seq_cds;
	print AXT_CDS $output_axt_cds;
    }	



sub convert_codon{
	my $codon = shift;
	my $ref_codon = shift;
	my $phase = shift;
	my @all;
	
	my $base = substr($codon,$phase,1);
	foreach my $replace (@{$Abbrev{$base}}) {
		my $new_codon = $ref_codon;
		substr($new_codon,$phase,1) = $replace;
		push @all,$new_codon;
	}

	return @all;
}



##caculate codon and phase
#############################################
sub codon_phase {
	my $pos = shift;

	my $phase = ($pos-1)%3;
	my $codon = ($pos-1-$phase)/3 + 1;
	
	##print "$codon,$phase,,\n";
	return ($codon,$phase);

}



sub different_base {
	my ($ref,$rice) = @_;
	my $base;
#print"$rice";exit;
	if ($rice =~ /[ACGT]/) {
		$base = $rice;
	}elsif($rice =~ /[MRWSYK]/){
		$base = $Abbrev{$rice}[0] if($Abbrev{$rice}[0] ne $ref);
		$base = $Abbrev{$rice}[1] if($Abbrev{$rice}[1] ne $ref);
		
	}else{
#print"$rice";exit;
	die "different_base unknown complex character, please check the snp data";
	}
	##print "$ref,$rice,$base;\n";
	return $base;
}

##get the status of a snp
#############################################
sub hom_het {
	my ($ref,$rice) = @_;
	
	return "same" if($ref eq $rice);
	
	my $status;
	

	if ($rice =~ /[ACGT]/) {
		$status = "Hom";
	}elsif($rice =~ /[MRWSYK]/){
		if ($Abbrev{$rice}[0] ne $ref && $Abbrev{$rice}[1] ne $ref) {
			$status = "Het-two";
		}else{
			$status = "Het-one";
		}
	}else{
		die "hom_het unknown complex character, please check the snp data";
	}

	return $status;
}

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


##complement and reverse the given sequence
#usage: Complement_Reverse(\$seq);
#############################################
sub Complement_Reverse{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
		$$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
		$$seq_p=reverse($$seq_p);  
	}
}
#############################################
