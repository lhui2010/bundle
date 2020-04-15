#!/usr/bin/perl -w

#usage: calc_RP.PL xx.bowtie xx.gff >result

#change dir back to file, gff is read into 24 hash tables
#2010/8/24

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
#1  2  3  4...  5...
#gene  counts  length  line_q  line_gff
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

$inittime=time;


if (@ARGV==0)
{
  print "usage:anno_bowtie.pl  *.bowtie XX.gff\n";
  exit;
}

@filenames=@ARGV;

$gff_name=pop @filenames;

&readgff;

$gfftime=time;
$duringtime=$gfftime-$inittime;

print "Begin...\nReading gff cost:$duringtime\n";

$lasttime=$gfftime;


#begin main

for$filename(@filenames)
{
open BOWTIE, $filename;

@e=split /\./,$filename;

$bsname=shift@e;

open OUT, ">$bsname".".anno_in_genome";

while(<BOWTIE>)
{
  chomp;

  @element_q=split /\s+/, $_;  #split bowtie
  push(@element_q, "")  if(@element_q == 8);  #no mismatch
  $line_q=join"\t", @element_q;  #storing the orignal line in .bowtie for print

  $reads_name=shift(@element_q);  #shift from end
  next unless($reads_name);  #for exception  

  $reads_number=shift(@element_q);
  $strand_q=shift(@element_q);#+
  $chr_q=shift(@element_q);#chrX

  if($strand_q eq "+")
  {
    $chr_q.="_1";
  }
  else
  {
    $chr_q.="_0";
  } 

  $pos_q=shift(@element_q);#132444
  $seq_q=shift(@element_q);
  $reads_len=length($seq_q);

#  close GFF;
#  open GFF, "$gffdir\/$chr_q$strand_q" ||warn"$!:$chr_q$strand_q\n";

  @rst=&withingff;
  if(@rst)#containing gene's name and RPB
  {
    for $rst(@rst)
    {
      print OUT "$rst\n";
    }
  }
  else
  {
    $strand_no = 0;
    $strand_no = 1 if($strand_q eq '+');
    print OUT "ncRNA.$chr_q".".$strand_no.$pos_q\t$reads_number\t$reads_len\t$line_q\n";
  }

}

close BOWTIE;
close OUT;

$thistime=time-$lasttime;
print"$bsname\n---------------------------------\n$thistime s\n\n";

$lasttime=time;

} 

$totaltime=time-$inittime;
print"\nTotal cost: $totaltime\n";

#---------------------------------------subs---------------------------------------------

#func:
#gff into hash table
#additional feature:identify 5utr intron exon and 3utr 
sub readgff
{
#chr1_1/0 means sense strand(1) of chr1 and anti..
#hash tables consists of chr1_1/0-chrM_1/0,and postion(xx-yy)is the key

  open GFF, "$gff_name" or die "can't open GFF file:$gff_name\n";

  my $lastname="";

  while(<GFF>)
  {
    chomp;
    @element_GFF=split;
    $gff_line=join"\t", @element_GFF;#storing the original line in gff for print()

    $chr=$element_GFF[0];
    next if(length($chr)>5);#ignore chr2389r03

    $chr="chrM" if($chr eq "chrMT");#consistent with genome

  #  $type=$element_GFF[2];

    $id=$element_GFF[-1];#-1
    @ids=split/;/,$id;
    shift(@ids);
    $id=shift(@ids);#name"XX"
    @ids=split'=',$id;#YY="XX"
    shift(@ids);
    $name=shift(@ids);#XX

    $start=$element_GFF[3];
    $end=$element_GFF[4];
      
    next if($name eq $lastname);

    $lastname=$name;

    if($element_GFF[6] eq "+")
    {
      $chr.="_1";#sense strand
    }
    else
    {
      $chr.="_0";
    }
    &buildhash();    
  }
  close GFF;  
}

sub buildhash
{
  for($chr)
  {
    if(/chr1_1/){$chr1_1{$start."-".$end}=$gff_line; push @tmp_chr1_1, $start." ".$end}
    elsif(/chr2_1/){$chr2_1{$start."-".$end}=$gff_line; push @tmp_chr2_1, $start." ".$end}
    elsif(/chr3_1/){$chr3_1{$start."-".$end}=$gff_line; push @tmp_chr3_1, $start." ".$end}
    elsif(/chr4_1/){$chr4_1{$start."-".$end}=$gff_line; push @tmp_chr4_1, $start." ".$end}
    elsif(/chr5_1/){$chr5_1{$start."-".$end}=$gff_line; push @tmp_chr5_1, $start." ".$end}
    elsif(/chr6_1/){$chr6_1{$start."-".$end}=$gff_line; push @tmp_chr6_1, $start." ".$end}
    elsif(/chr7_1/){$chr7_1{$start."-".$end}=$gff_line; push @tmp_chr7_1, $start." ".$end}
    elsif(/chr8_1/){$chr8_1{$start."-".$end}=$gff_line; push @tmp_chr8_1, $start." ".$end}
    elsif(/chr9_1/){$chr9_1{$start."-".$end}=$gff_line; push @tmp_chr9_1, $start." ".$end}
    elsif(/chr10_1/){$chr10_1{$start."-".$end}=$gff_line; push @tmp_chr10_1, $start." ".$end}
    elsif(/chr11_1/){$chr11_1{$start."-".$end}=$gff_line; push @tmp_chr11_1, $start." ".$end}
    elsif(/chr12_1/){$chr12_1{$start."-".$end}=$gff_line; push @tmp_chr12_1, $start." ".$end}
    elsif(/chr13_1/){$chr13_1{$start."-".$end}=$gff_line; push @tmp_chr13_1, $start." ".$end}
    elsif(/chr14_1/){$chr14_1{$start."-".$end}=$gff_line; push @tmp_chr14_1, $start." ".$end}
    elsif(/chr15_1/){$chr15_1{$start."-".$end}=$gff_line; push @tmp_chr15_1, $start." ".$end}
    elsif(/chr16_1/){$chr16_1{$start."-".$end}=$gff_line; push @tmp_chr16_1, $start." ".$end}
    elsif(/chr17_1/){$chr17_1{$start."-".$end}=$gff_line; push @tmp_chr17_1, $start." ".$end}
    elsif(/chr18_1/){$chr18_1{$start."-".$end}=$gff_line; push @tmp_chr18_1, $start." ".$end}
    elsif(/chr19_1/){$chr19_1{$start."-".$end}=$gff_line; push @tmp_chr19_1, $start." ".$end}
    elsif(/chr20_1/){$chr20_1{$start."-".$end}=$gff_line; push @tmp_chr20_1, $start." ".$end}
    elsif(/chr21_1/){$chr21_1{$start."-".$end}=$gff_line; push @tmp_chr21_1, $start." ".$end}
    elsif(/chr22_1/){$chr22_1{$start."-".$end}=$gff_line; push @tmp_chr22_1, $start." ".$end}
    elsif(/chrX_1/){$chrX_1{$start."-".$end}=$gff_line; push @tmp_chrX_1, $start." ".$end}
    elsif(/chrY_1/){$chrY_1{$start."-".$end}=$gff_line; push @tmp_chrY_1, $start." ".$end}
    elsif(/chrM_1/){$chrM_1{$start."-".$end}=$gff_line; push @tmp_chrM_1, $start." ".$end}
    elsif(/chr1_0/){$chr1_0{$start."-".$end}=$gff_line; push @tmp_chr1_0, $start." ".$end}
    elsif(/chr2_0/){$chr2_0{$start."-".$end}=$gff_line; push @tmp_chr2_0, $start." ".$end}
    elsif(/chr3_0/){$chr3_0{$start."-".$end}=$gff_line; push @tmp_chr3_0, $start." ".$end}
    elsif(/chr4_0/){$chr4_0{$start."-".$end}=$gff_line; push @tmp_chr4_0, $start." ".$end}
    elsif(/chr5_0/){$chr5_0{$start."-".$end}=$gff_line; push @tmp_chr5_0, $start." ".$end}
    elsif(/chr6_0/){$chr6_0{$start."-".$end}=$gff_line; push @tmp_chr6_0, $start." ".$end}
    elsif(/chr7_0/){$chr7_0{$start."-".$end}=$gff_line; push @tmp_chr7_0, $start." ".$end}
    elsif(/chr8_0/){$chr8_0{$start."-".$end}=$gff_line; push @tmp_chr8_0, $start." ".$end}
    elsif(/chr9_0/){$chr9_0{$start."-".$end}=$gff_line; push @tmp_chr9_0, $start." ".$end}
    elsif(/chr10_0/){$chr10_0{$start."-".$end}=$gff_line; push @tmp_chr10_0, $start." ".$end}
    elsif(/chr11_0/){$chr11_0{$start."-".$end}=$gff_line; push @tmp_chr11_0, $start." ".$end}
    elsif(/chr12_0/){$chr12_0{$start."-".$end}=$gff_line; push @tmp_chr12_0, $start." ".$end}
    elsif(/chr13_0/){$chr13_0{$start."-".$end}=$gff_line; push @tmp_chr13_0, $start." ".$end}
    elsif(/chr14_0/){$chr14_0{$start."-".$end}=$gff_line; push @tmp_chr14_0, $start." ".$end}
    elsif(/chr15_0/){$chr15_0{$start."-".$end}=$gff_line; push @tmp_chr15_0, $start." ".$end}
    elsif(/chr16_0/){$chr16_0{$start."-".$end}=$gff_line; push @tmp_chr16_0, $start." ".$end}
    elsif(/chr17_0/){$chr17_0{$start."-".$end}=$gff_line; push @tmp_chr17_0, $start." ".$end}
    elsif(/chr18_0/){$chr18_0{$start."-".$end}=$gff_line; push @tmp_chr18_0, $start." ".$end}
    elsif(/chr19_0/){$chr19_0{$start."-".$end}=$gff_line; push @tmp_chr19_0, $start." ".$end}
    elsif(/chr20_0/){$chr20_0{$start."-".$end}=$gff_line; push @tmp_chr20_0, $start." ".$end}
    elsif(/chr21_0/){$chr21_0{$start."-".$end}=$gff_line; push @tmp_chr21_0, $start." ".$end}
    elsif(/chr22_0/){$chr22_0{$start."-".$end}=$gff_line; push @tmp_chr22_0, $start." ".$end}
    elsif(/chrX_0/){$chrX_0{$start."-".$end}=$gff_line; push @tmp_chrX_0, $start." ".$end}
    elsif(/chrY_0/){$chrY_0{$start."-".$end}=$gff_line; push @tmp_chrY_0, $start." ".$end}
    elsif(/chrM_0/){$chrM_0{$start."-".$end}=$gff_line; push @tmp_chrM_0, $start." ".$end}
    else{}
  }

#sort
  @chr1_1=sort by_start@tmp_chr1_1;
  @chr2_1=sort by_start@tmp_chr2_1;
  @chr3_1=sort by_start@tmp_chr3_1;
  @chr4_1=sort by_start@tmp_chr4_1;
  @chr5_1=sort by_start@tmp_chr5_1;
  @chr6_1=sort by_start@tmp_chr6_1;
  @chr7_1=sort by_start@tmp_chr7_1;
  @chr8_1=sort by_start@tmp_chr8_1;
  @chr9_1=sort by_start@tmp_chr9_1;
  @chr10_1=sort by_start@tmp_chr10_1;
  @chr11_1=sort by_start@tmp_chr11_1;
  @chr12_1=sort by_start@tmp_chr12_1;
  @chr13_1=sort by_start@tmp_chr13_1;
  @chr14_1=sort by_start@tmp_chr14_1;
  @chr15_1=sort by_start@tmp_chr15_1;
  @chr16_1=sort by_start@tmp_chr16_1;
  @chr17_1=sort by_start@tmp_chr17_1;
  @chr18_1=sort by_start@tmp_chr18_1;
  @chr19_1=sort by_start@tmp_chr19_1;
  @chr20_1=sort by_start@tmp_chr20_1;
  @chr21_1=sort by_start@tmp_chr21_1;
  @chr22_1=sort by_start@tmp_chr22_1;
  @chrX_1=sort by_start@tmp_chrX_1;
  @chrY_1=sort by_start@tmp_chrY_1;
  @chrM_1=sort by_start@tmp_chrM_1;
  @chr1_0=sort by_start@tmp_chr1_0;
  @chr2_0=sort by_start@tmp_chr2_0;
  @chr3_0=sort by_start@tmp_chr3_0;
  @chr4_0=sort by_start@tmp_chr4_0;
  @chr5_0=sort by_start@tmp_chr5_0;
  @chr6_0=sort by_start@tmp_chr6_0;
  @chr7_0=sort by_start@tmp_chr7_0;
  @chr8_0=sort by_start@tmp_chr8_0;
  @chr9_0=sort by_start@tmp_chr9_0;
  @chr10_0=sort by_start@tmp_chr10_0;
  @chr11_0=sort by_start@tmp_chr11_0;
  @chr12_0=sort by_start@tmp_chr12_0;
  @chr13_0=sort by_start@tmp_chr13_0;
  @chr14_0=sort by_start@tmp_chr14_0;
  @chr15_0=sort by_start@tmp_chr15_0;
  @chr16_0=sort by_start@tmp_chr16_0;
  @chr17_0=sort by_start@tmp_chr17_0;
  @chr18_0=sort by_start@tmp_chr18_0;
  @chr19_0=sort by_start@tmp_chr19_0;
  @chr20_0=sort by_start@tmp_chr20_0;
  @chr21_0=sort by_start@tmp_chr21_0;
  @chr22_0=sort by_start@tmp_chr22_0;
  @chrX_0=sort by_start@tmp_chrX_0;
  @chrY_0=sort by_start@tmp_chrY_0;
  @chrM_0=sort by_start@tmp_chrM_0;

}

sub by_start
{
@a=split/\s/, $a;
@b=split/\s/, $b;
$a[0]<=>$b[0];
}



 
#following global var was used: 
#$reads_number
#$chr_q $strand_q $pos_q $reads_len
#return gene's name and RPB
sub withingff
{
  my @rstl;
  my $lastname="";#incase overlapping gene appear
  #select
  for($chr_q)
  {
    if(/chr1_1/){$phash=\%chr1_1;$parray=\@chr1_1;}
    elsif(/chr2_1/){$phash=\%chr2_1;$parray=\@chr2_1;}
    elsif(/chr3_1/){$phash=\%chr3_1;$parray=\@chr3_1;}
    elsif(/chr4_1/){$phash=\%chr4_1;$parray=\@chr4_1;}
    elsif(/chr5_1/){$phash=\%chr5_1;$parray=\@chr5_1;}
    elsif(/chr6_1/){$phash=\%chr6_1;$parray=\@chr6_1;}
    elsif(/chr7_1/){$phash=\%chr7_1;$parray=\@chr7_1;}
    elsif(/chr8_1/){$phash=\%chr8_1;$parray=\@chr8_1;}
    elsif(/chr9_1/){$phash=\%chr9_1;$parray=\@chr9_1;}
    elsif(/chr10_1/){$phash=\%chr10_1;$parray=\@chr10_1;}
    elsif(/chr11_1/){$phash=\%chr11_1;$parray=\@chr11_1;}
    elsif(/chr12_1/){$phash=\%chr12_1;$parray=\@chr12_1;}
    elsif(/chr13_1/){$phash=\%chr13_1;$parray=\@chr13_1;}
    elsif(/chr14_1/){$phash=\%chr14_1;$parray=\@chr14_1;}
    elsif(/chr15_1/){$phash=\%chr15_1;$parray=\@chr15_1;}
    elsif(/chr16_1/){$phash=\%chr16_1;$parray=\@chr16_1;}
    elsif(/chr17_1/){$phash=\%chr17_1;$parray=\@chr17_1;}
    elsif(/chr18_1/){$phash=\%chr18_1;$parray=\@chr18_1;}
    elsif(/chr19_1/){$phash=\%chr19_1;$parray=\@chr19_1;}
    elsif(/chr20_1/){$phash=\%chr20_1;$parray=\@chr20_1;}
    elsif(/chr21_1/){$phash=\%chr21_1;$parray=\@chr21_1;}
    elsif(/chr22_1/){$phash=\%chr22_1;$parray=\@chr22_1;}
    elsif(/chrX_1/){$phash=\%chrX_1;$parray=\@chrX_1;}
    elsif(/chrY_1/){$phash=\%chrY_1;$parray=\@chrY_1;}
    elsif(/chrM_1/){$phash=\%chrM_1;$parray=\@chrM_1;}
    elsif(/chr1_0/){$phash=\%chr1_0;$parray=\@chr1_0;}
    elsif(/chr2_0/){$phash=\%chr2_0;$parray=\@chr2_0;}
    elsif(/chr3_0/){$phash=\%chr3_0;$parray=\@chr3_0;}
    elsif(/chr4_0/){$phash=\%chr4_0;$parray=\@chr4_0;}
    elsif(/chr5_0/){$phash=\%chr5_0;$parray=\@chr5_0;}
    elsif(/chr6_0/){$phash=\%chr6_0;$parray=\@chr6_0;}
    elsif(/chr7_0/){$phash=\%chr7_0;$parray=\@chr7_0;}
    elsif(/chr8_0/){$phash=\%chr8_0;$parray=\@chr8_0;}
    elsif(/chr9_0/){$phash=\%chr9_0;$parray=\@chr9_0;}
    elsif(/chr10_0/){$phash=\%chr10_0;$parray=\@chr10_0;}
    elsif(/chr11_0/){$phash=\%chr11_0;$parray=\@chr11_0;}
    elsif(/chr12_0/){$phash=\%chr12_0;$parray=\@chr12_0;}
    elsif(/chr13_0/){$phash=\%chr13_0;$parray=\@chr13_0;}
    elsif(/chr14_0/){$phash=\%chr14_0;$parray=\@chr14_0;}
    elsif(/chr15_0/){$phash=\%chr15_0;$parray=\@chr15_0;}
    elsif(/chr16_0/){$phash=\%chr16_0;$parray=\@chr16_0;}
    elsif(/chr17_0/){$phash=\%chr17_0;$parray=\@chr17_0;}
    elsif(/chr18_0/){$phash=\%chr18_0;$parray=\@chr18_0;}
    elsif(/chr19_0/){$phash=\%chr19_0;$parray=\@chr19_0;}
    elsif(/chr20_0/){$phash=\%chr20_0;$parray=\@chr20_0;}
    elsif(/chr21_0/){$phash=\%chr21_0;$parray=\@chr21_0;}
    elsif(/chr22_0/){$phash=\%chr22_0;$parray=\@chr22_0;}
    elsif(/chrX_0/){$phash=\%chrX_0;$parray=\@chrX_0;}
    elsif(/chrY_0/){$phash=\%chrY_0;$parray=\@chrY_0;}
    elsif(/chrM_0/){$phash=\%chrM_0;$parray=\@chrM_0;}
  }

  
  for($i=0; $i<@$parray; $i++)
  {
    $_ =$$parray[$i];

    @start_end=split;

    $start=shift @start_end;
    $end=shift @start_end;
    
    if($pos_q>= $start && $pos_q<=$end)
    {
      $gff_line=$$phash{$start."-".$end};

      @element_GFF=split/\t/, $gff_line;
      
      $id=pop @element_GFF;#-1
      @ids=split/;/,$id;
      shift(@ids);
      $id=shift(@ids);#name"XX"
      @ids=split'=',$id;#YY="XX"
      shift(@ids);
      $name=shift(@ids);#XX

      next if($name eq $lastname);
      push @rstl,"$name\t$reads_number\t$reads_len\t$line_q\t$gff_line";
      $lastname=$name;

    }
    #assuming genes in gff are sorted by position
    if($pos_q<$start)
    {
      return @rstl;
    }
  }
  @rstl;
}
