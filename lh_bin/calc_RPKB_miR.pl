#!/usr/bin/perl -w

#calc_RPB will calculate RPB(readsnum * readslen/ncRNAlen) from bowtie default output file
# and combining GFF annotation to the map result.
#Thus two new column was added ab ovo: Gene's name, RPB, and reads length
#The sitution that multiple reads mapped to the same gens was ignored in 
#this version, so there may be duplicated genes in the result

#usage: calc_RPB XX.gff XX(bowtie's output) -n/a >result
#n for novel output
#a for annotated output
#Input: XX.gff, XX.bowtie
#Output: print result, XX.bowtie_in_XX.gff(number of reads in .bowtie, number of matched reads in .gff)

print "usage:calc_RPB XX.gff XX.bowtie\n" if (@ARGV==0);

@filenames=@ARGV;

$option=$ARGV[2];

open BOWTIE, $filenames[1];

while(<BOWTIE>)
{
  chomp;
  $line_q=$_;

  $number_reads++;

  open GFF, $filenames[0];

  @element_q=split /\s+/, $_;#split bowtie

#  $reads_name=
  shift(@element_q);#shift from end
  $reads_number=shift(@element_q);
  $strand_q=shift(@element_q);#+
  $chr_q=shift(@element_q);#chrX
  $pos_q=shift(@element_q);#132444
  $seq_q=shift(@element_q);
  $reads_len=length($seq_q);

  if($option eq '-n')
  {
	$novel++;
	print"ncRNA.$novel\t$reads_number\t$reads_len\t$line_q\n";
  }
  elsif($option eq '-a')#containing gene's name and RPB
  {
    $rst=&withingff;
    if($rst)
    {
      print "$rst\t$reads_len\t$line_q\n";
      $number_matched++;
    }
  }

  close GFF;
}

open STAT, ">".$ARGV[1].".in_".$ARGV[0];
print STAT "number_matched:$number_matched\nnumber_reads:$number_reads\n";
  
#following global var was used: 
#$reads_number
#$chr_q $strand_q $pos_q $reads_len
#return gene's name and RPB
sub withingff
{
	while(<GFF>)
	{
		chomp;
		@element_GFF=split/\s+/,$_;
		$chr=shift(@element_GFF);#1
		$chr="chr".$chr;
		if($chr_q eq $chr)
		{
			shift(@element_GFF);#2
			$type=shift(@element_GFF);#3
			if($type eq "exon" || $type eq "miRNA")
			{
				$start=shift(@element_GFF);#4
				$end=shift(@element_GFF);#5
				if($pos_q>= $start && $pos_q<=$end)
				{
					shift(@element_GFF);#6
					$strand=shift(@element_GFF);#7
					if($strand_q eq $strand)#found finally
					{
						shift(@element_GFF);#8
						shift(@element_GFF);#9
						$id=shift(@element_GFF);#10
#print"\nDEBUG:$id\n";
#r, why miRb's gff is diff from ncbi's
#						@ids=split/;/,$id;
#						shift(@ids);
#						$id=shift(@ids);#name"XX"

						@ids=split'=',$id;#YY="XX"
						shift(@ids);
						$id=shift(@ids);#"XX"

						@ids=split'\"',$id;
						$name=$ids[1];

						$RPB=int($reads_number*$reads_len/($end-$start));

						return "$name\t$reads_number";
					}
				}
			}
		}
	}
	return 0;
}
