#!/usr/bin/perl-w
use strict;
use warnings;
use FileHandle;
#modified by liuhui
#Wed Apr 13 19:38:33 CST 2011

#explanation:this program is edited to
#edit by Lvjun;   Tue Apr 27 17:03:19 CST 2010
#Version 1.0    hewm@genomics.org.cn

die  "Version 1.0\t2010-04-27;\nUsage: $0 <InPut_1:CNS_filelist><OutDir>\n" unless (@ARGV ==2);

#############Befor  Start  ,open the files #########################

open     INFilist,"$ARGV[0]"  or die "input file can't open $!" ;
open     OutFile,">$ARGV[1]" or die "output file can't open $!" ;


############ Do what you want to do #####################

my $reads_count_threshold = 3;

my @file = ();
my %fh = ();
my $file_count = 0;

        while(<INFilist>)
        {
                chomp ;
		$file_count++;
		$file[$file_count] = $_;
		$fh{$file_count} = FileHandle->new("gzip -dc $_|");#new("gzip -dc $_|");
        }
	
close INFilist;
	my %het_gen = ("AC" => "M","CA" => "M","AT" => "W", "TA" => "W", "AG" => "R", "GA" => "R", "CT" => "Y", "TC" => "Y", "CG" => "S" , "GC" => "S", "GT" => "K", "TG" => "K", "AA" => "A", "TT" => "T", "CC" => "C", "GG" => "G");
	my %het_hom = ("A" => 1, "T" =>1, "G" => 1, "C" => 1, "M" => 0, "W" => 0, "R" => 0, "Y" => 0, "S" => 0, "K" => 0);
	while($_=$fh{1}->getline())
#	while(<>)
	{
		chomp;
	#	print "$_\n";exit;
		my @inf = split;
		my $gen;

		if ($inf[13] == 0) {$inf[3] = "-"; print OutFile "$inf[0]\t$inf[1]\t$inf[2]\t$inf[3] " ;}
	#	print OutFile "$inf[0]\t$inf[1]$inf[3]\t$inf[4]\t$inf[7]+$inf[11]\t$inf[13]";

#heterozygous site
		elsif( $inf[6] >= 20 && $inf[10] >= 20 && $inf[8] >= $reads_count_threshold && $inf[12] >=$reads_count_threshold)  ### assign the hetezygous genotype
		{
			$gen = "$inf[5]"."$inf[9]";
 			print OutFile "$inf[0]\t$inf[1]\t$inf[2]\t$het_gen{$gen} " ; 
		}

#not hetrozygous, best is best
		elsif( $inf[6] >=20 && $inf[8] >= $reads_count_threshold && $inf[10] < 20 && $het_hom{$inf[3]} ==1)		
		{	print OutFile "$inf[0]\t$inf[1]\t$inf[2]\t$inf[5] ";}	### assign the homozygous genotype

#secondbest is best
		elsif( $inf[6] < 20 && $inf[10] >= 20 && $inf[12] >=$reads_count_threshold && $het_hom{$inf[3]} ==1)
		{	print OutFile "$inf[0]\t$inf[1]\t$inf[2]\t$inf[9] ";} ### assign the homozygous genotype

#not homo, print hetero-one
		elsif( $inf[4] >= 20 && $inf[8]+$inf[12] >= $reads_count_threshold)		### if one base quality larger than 20 but the consensus is hetezygous and quality still more than 20 then it's hetezygous, or both the two base quality less than 20 but the consensus genotype larger than 20 then it's hetezygous
		{	print OutFile "$inf[0]\t$inf[1]\t$inf[2]\t$inf[3] ";}
		else
		{
			$inf[3] = "-"; print OutFile "$inf[0]\t$inf[1]\t$inf[2]\t$inf[3] ";
		}

		for(my $i = 2; $i <= $file_count; $i++)
		{
			my $rec = $fh{$i}->getline();
			chomp $rec;
			#print "$rec\n";
			my @temp = split /\t/,$rec;
			#print "@temp\n"; exit;
			if($inf[1] != $temp[1])
				{
					die "bad at the file $i pos $temp[1] and 1 pos $inf[1]\n";
				}
			my $gen;
               		if ($temp[13] == 0) {$temp[3] = "-"; print OutFile "$temp[3] " ;}
	        #       print OutFile "$temp[0]\t$temp[1]$temp[3]\t$temp[4]\t$temp[7]+$temp[11]\t$temp[13]";
        	        elsif( $temp[6] >= 20 && $temp[10] >= 20 && $temp[8] >= $reads_count_threshold && $temp[12] >=$reads_count_threshold)  ### assign the hetezygous genotype
        	        {
        	        $gen = "$temp[5]"."$temp[9]";
        	        print OutFile "$het_gen{$gen} ";
        	        }
        	        elsif( $temp[6] >=20 && $temp[8] >= $reads_count_threshold && $temp[10] < 20 && $het_hom{$temp[3]} ==1)
               		 {       print OutFile "$temp[5] ";} ### assign the homozygous genotype
               		 elsif( $temp[6] < 20 && $temp[10] >= 20 && $temp[12] >=$reads_count_threshold && $het_hom{$temp[3]} ==1)
               		 {       print OutFile "$temp[9] ";} ### assign the homozygous genotype
               		 elsif( $temp[4] >= 20 && $temp[8]+$temp[12] >= $reads_count_threshold)           ### if one base quality larger than 20 but the consensus is hetezygous and quality still more than 20 then it's hetezygous, or both the two base quality less than 20 but the consensus genotype larger than 20 then it's hetezygous
               		 {   
				#$gen = "$temp[5]"."$temp[9]";
			     	print OutFile "$temp[3] ";
			}
             		   else
                	{$temp[3] = "-"; print OutFile "$temp[3] ";}
		}
		print OutFile "\n";
	}
	
close OutFile ;

