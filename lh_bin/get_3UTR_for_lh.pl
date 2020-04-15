#!/usr/bin/perl-w
use strict;
use warnings;
#explanation:this program is edited to
#edit by Lvjun;   Tue Apr 27 17:03:19 CST 2010
#Version 1.0    hewm@genomics.org.cn

die  "Version 1.0\t2010-04-27;\nUsage: $0 <InPut_1><InPut_2><OutDir>\n" unless (@ARGV ==3);

#############Befor  Start  ,open the files #########################

open     INFile,"$ARGV[0]"  || die "input file can't open $!" ;
open     INFile2,"$ARGV[1]"  || die "input file can't open $!" ;
open     OutFile,">$ARGV[2]" || die "output file can't open $!" ;

my %UTR3_gene_info;

############ Do what you want to do #####################

        while(<INFile>)
        {
                chomp ;
                my @inf=split ;
		my $info = join "\t",@inf[0..$#inf];
		for(my $i = $inf[2]; $i<=$inf[3]; $i++)
		{
			if(!exists $UTR3_gene_info{$i})
			{ $UTR3_gene_info{$i} = $info; }
		}


        }
	while(<INFile2>)
	{
		chomp;
		my $temp = $_;
		if(exists $UTR3_gene_info{$temp})
		{ print OutFile "$temp\t$UTR3_gene_info{$temp}\n";}
	
	}

close INFile;
close OutFile ;

