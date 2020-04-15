#!/usr/bin/perl-w
use strict;
use warnings;

#modified by liuhui to adjust to .gz input
#Sun Apr 17 12:40:29 CST 2011

#explanation:this program is edited to
#edit by Lvjun;   Tue Apr 27 17:03:19 CST 2010
#Version 1.0    hewm@genomics.org.cn

die  "Version 1.0\t2010-04-27;\nUsage: $0 <InPut_1:start_codon_region_file><InPut_2:RSAs_file><OutDir:RSA_in_start_codon>\n" unless (@ARGV ==3);

#############Befor  Start  ,open the files #########################

open     START,"$ARGV[0]"  || die "input file can't open $!" ;
open     RSA,"gzip -dc $ARGV[1]|"  || die "input file can't open $!" ;
open     OutFile,">$ARGV[2]" || die "output file can't open $!" ;


############ Do what you want to do #####################
my %inter;
        while(<START>)
        {
                chomp ;
                my @inf=split /\t/;
		for(my $i = $inf[2]; $i <= $inf[3]; $i++)	
		{
			$inter{$i} = 1; 
		}
        }
	while(<RSA>)
	{
		chomp;
		my @temp = split /\t/;
		next if(!exists $inter{$temp[1]});
		my $t = join "\t",@temp[0..$#temp];
		print OutFile "$t\n";
	
	}

#close INTER;
close RSA;
close OutFile ;

