#!/usr/bin/perl -w

#usage: calc_RPB gff_directory XX(bowtie's output) >result

#change expr, calc_RPB gff_directory xx.bowtie >result
#change name from calc_RPKB.pl to anno_bowtie.pl
#try to fix "open on closed file handler" bug
#21:04 2010/8/17

#error again, gff seems never been opened, kinds of ncRNA...
#masking all close gff now
#16:12 2010/8/4

#fixed some logical error, append "" to input array if it's
#perfect match
#14:37 1010/8/4

#new change: appending matching line in .gff to the end"
#format specifation:
#1	2	3	4...	5...
#gene	counts	length	line_q	line_gff
#
#gene:gene's name
#counts: reads counts of that gene
#length: length of querying reads
#line_q: the rest line of querying bowtie file(in default format)
#line_gff: total line from line_gff file(in gff3 format)
#
#ps: no more shift in sub withingff()
#13:42 2010/8/1
#pps: leave shift alone, copy to a new array instead
#13:45 2010/8/1

#calc_RPB will calculate RPB(readsnum * readslen/ncRNAlen) from bowtie default output file
# and combining GFF annotation to the map result.
#Thus two new column was added ab ovo: Gene's name, RPB, and reads length
#The sitution that multiple reads mapped to the same gens was ignored in 
#this version, so there may be duplicated genes in the result

#new: when supplies XX.gff, calc_RPKB.pl will actually go for _XX.gff dir to search for annotation

if (@ARGV==0)
{
  print "usage:calc_RPB XX.gff XX.bowtie\n";
  exit;
}

@filenames=@ARGV;


open BOWTIE, $filenames[1];

$gffdir=$filenames[0];

while(<BOWTIE>)
{
  chomp;

  @element_q=split /\s+/, $_;	#split bowtie
  push(@element_q, "")  if(@element_q == 8);	#no mismatch
  $line_q=join"\t", @element_q;	#storing the orignal line in .bowtie for print

  $reads_name=shift(@element_q);	#shift from end
  next unless($reads_name);	#for exception  

  $reads_number=shift(@element_q);
  $strand_q=shift(@element_q);#+
  $chr_q=shift(@element_q);#chrX
  $pos_q=shift(@element_q);#132444
  $seq_q=shift(@element_q);
  $reads_len=length($seq_q);

#  close GFF;
#  open GFF, "$gffdir\/$chr_q$strand_q" ||warn"$!:$chr_q%strand_q\n";

  $rst=&withingff;
  if($rst)#containing gene's name and RPB
	{
	  print "$rst\t$reads_len\t$line_q\t$gff_print\n";
	}
	else
	{
		$strand_no = 0;
		$strand_no = 1 if($strand_q eq '+');
		print"ncRNA.$chr_q".".$strand_no.$pos_q\t$reads_number\t$reads_len\t$line_q\n";
	}
#  close GFF;
}

close BOWTIE;
  
#following global var was used: 
#$reads_number
#$chr_q $strand_q $pos_q $reads_len
#return gene's name and RPB
sub withingff
{
	open GFF, "$gffdir\/$chr_q$strand_q" or die"can't open $gffdir\/$chr_q$strand_q\n$!\n";
#return 0;

	while(<GFF>)
	{
		chomp;

		@element_GFF=split;
		$gff_print=join"\t", @element_GFF;#storing the original line in gff for print()

		shift(@element_GFF);#1
		shift(@element_GFF);#2
		shift(@element_GFF);#3

		$start=shift(@element_GFF);#4
		$end=shift(@element_GFF);#5
		if($pos_q>= $start && $pos_q<=$end)
		{
			$id=pop @element_GFF;#-1
			@ids=split/;/,$id;
			shift(@ids);
			$id=shift(@ids);#name"XX"
			@ids=split'=',$id;#YY="XX"
			shift(@ids);
			$name=shift(@ids);#XX

			close GFF;
			return "$name\t$reads_number";

		}
	}
	close GFF;
	0;
}
