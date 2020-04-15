#!/usr/bin/perl -w
use strict;


#TODO: current version only support single fasta file
#So we need to update $seq to %hash{chr1} to support multi-fasta file


if(@ARGV < 4)
{
	print "function:\n both start nt and end nt will be output\nusage:\n perl -w $0 start end chromosome_file snp_file\n";
	exit;
}

my $start0 = shift @ARGV;
my $end0 = shift @ARGV;
my $chr_file = shift;
my $snp_file = shift;


####################################################
#  ____                _   _____         _         #
# |  _ \ ___  __ _  __| | |  ___|_ _ ___| |_ __ _  #
# | |_) / _ \/ _` |/ _` | | |_ / _` / __| __/ _` | #
# |  _ <  __/ (_| | (_| | |  _| (_| \__ \ || (_| | #
# |_| \_\___|\__,_|\__,_| |_|  \__,_|___/\__\__,_| #
#                                                  #
####################################################

#will read whole fasta to $seq; then use substr($seq, 1,1) = xx to manipulate.

$/='>';
my $seq = "";

open FA, $chr_file or die;
while(<FA>)
{
	chomp;
	next if ($_ eq "");
	my $info;

	($info, $seq) = split "\n", $_, 2;
	$seq =~ s/\n//g;
#	$hash{$info} = $seq;
}

$/="\n";

($seq) = (">".$seq);





###############################################
#  ____                _   ____  _   _ ____   #
# |  _ \ ___  __ _  __| | / ___|| \ | |  _ \  #
# | |_) / _ \/ _` |/ _` | \___ \|  \| | |_) | #
# |  _ <  __/ (_| | (_| |  ___) | |\  |  __/  #
# |_| \_\___|\__,_|\__,_| |____/|_| \_|_|     #
#                                             #		
###############################################

#will read SNP specified segment; change ref in $seq to $snp;





#####################################################
#65 list
#@acc_list = qw/12883_AUS 45975_AUS 6307_AUS 8555_AUS 27762_IND 30416_IND 43545_IND 51250_IND 51300_IND 8231_IND 9148_IND 9177_IND 1107_TEJ 2540_TEJ 27630_TEJ 32399_TEJ 418_TEJ 55471_TEJ 8191_TEJ NP_TEJ 11010_TRJ 17757_TRJ 328_TRJ 38698_TRJ 43325_TRJ 43397_TRJ 43675_TRJ 50448_TRJ 66756_TRJ 8244_TRJ 12793_ARO 38994_ARO 9060_ARO 9062_ARO RA4952_ARO nivara_103407 nivara_105327 nivara_105705 nivara_105784 nivara_105879 nivara_106105 nivara_106154 nivara_106345 nivara_80470 nivara_89215 rufipogon_105426 rufipogon_105912 rufipogon_105958 rufipogon_105960 rufipogon_106161 rufipogon_106505 rufipogon_80506 rufipogon_81982 rufipogon_81991 rufipogon_Dongxiang rufipogon_Nepal rufipogon_P25 rufipogon_P46 rufipogon_P61 rufipogon_YJ/;


#####################################################
#50 list
my @acc_list = qw/12883_AUS 45975_AUS 6307_AUS 8555_AUS 27762_IND 26872_TRJ 30416_IND 43545_IND 51250_IND 51300_IND 8231_IND 9148_IND 9177_IND 1107_TEJ 27630_TEJ 32399_TEJ 2540_TEJ 418_TEJ 55471_TEJ 8191_TEJ NP_TEJ 11010_TRJ 17757_TRJ 25901_IND 328_TRJ 38698_TRJ 43325_TRJ 43397_IND 43675_TRJ 50448_TRJ 66756_TRJ 8244_TRJ 12793_ARO 38994_ARO 9060_ARO 9062_ARO RA4952_ARO 31856_V 60542_IV 6513_III rufipogon_105958 rufipogon_105960 rufipogon_Nepal rufipogon_P46 rufipogon_YJ nivara_105327 nivara_106105 nivara_106154 nivara_80470 nivara_89215/;

my $column_count = int(@acc_list);



###################################################
#indicate whether to reverse complement
#reverse will take place in last step
my $rc_flag = 0;


#$start0 and $end0 is original input, which also indicate strand direction
#while $start always < $end.
my ($start, $end) = ($start0, $end0);

if($start0 > $end0)
{
	$rc_flag =1;
	($start, $end) = ($end, $start);
}
####################################################






my %acc_fa; #final fasta store
my $ref_fa = substr($seq, $start, $end - $start +1); #$seq is 1-coordinate based

for  (my $i=0; $i<$column_count; $i++)
{
	$acc_fa{$acc_list[$i]} = $ref_fa;
}

####Main########################################################
#	format of input <>
#	chromosome07    435     G G G G G G G G G G G G G G G G A G G G G R A G G A A G A G G A G A A A A A A A G G G R G G G G G G G
#################################################################
open SNP, $snp_file or die;
my $chr;

while(<SNP>)
{
	s/-/N/g;
	chomp;
	my @e=split;

	$chr = shift @e;
	my $loci = shift @e;
	my $ref = shift @e;

	my $rel_loci = $loci - $start;

	next if($loci < $start);
	last if($loci > $end);

	for (my $i=0; $i<$column_count; $i++)
	{
		substr($acc_fa{$acc_list[$i]}, $rel_loci, 1) = $e[$i];
	}
}


for (my $i=0; $i<$column_count; $i++)
{
	print ">$acc_list[$i]\t$chr:$start0-$end0\n";
	&rc(\$acc_fa{$acc_list[$i]}) if($rc_flag);
	&formated_print (\$acc_fa{$acc_list[$i]});
}



###########################################################
#   ____        _     ____             _   _              #
#  / ___| _   _| |__ |  _ \ ___  _   _| |_(_)_ __   ___   #
#  \___ \| | | | '_ \| |_) / _ \| | | | __| | '_ \ / _ \  #
#   ___) | |_| | |_) |  _ < (_) | |_| | |_| | | | |  __/s #
#  |____/ \__,_|_.__/|_| \_\___/ \__,_|\__|_|_| |_|\___|  #
#                                                         #
###########################################################




sub rc
{
            my $seq_p=shift;
                if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
                                $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
                                        $$seq_p=reverse($$seq_p);
                                        }
}

sub formated_print
{
	        my $length = 70;
		        my $seq_p =shift;
			        if (ref($seq_p) eq 'SCALAR') {
					                for ( my $pos = 0 ; $pos < length($$seq_p) ; $pos += $length ) {
								                        print substr($$seq_p, $pos, $length), "\n";
											                }
													        }
}

