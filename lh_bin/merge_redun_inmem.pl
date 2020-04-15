#!/usr/bin/perl -w

#forget about #3#4, it can't be done, I'll use column #1 only instead
#2010/8/2

#diff: read all data to memory and compare within it
#1. read all #1 column to memory
#2. sort in memory, repeated one is assigned to a value in the tag array. count unique line
#3. split the file into two file(redun, non-redun)s according to the tag array
#4. compress non-redundant array

#function: merged repeated lines by the first column
#
#usage: merge_redun.pl tobemergedfile(can be multi) >cout_file
#
#input format:
#1.tobemergedfile
#11 xx	Risc
#22 zz	fds
#11 yy	Bisc
#
#output specification:
#1.tobemergedfile_noredun
#11 xx	Risc
#	 yy	Bisc
#22 zz	fds
#2.screen output
#filename
#count of no redun line
#count of total line

#detailed function:
#merge redundance.in_miR(splited by \s+) by name
#repeated element(#1)'s line will be regarded as once
#and stored in one array for indentifer and another
#Each unique name is printed only once;
# with its rest lines for output
#eg:
#arrayiden:(11,22,33)
#arrayprint:(11\txx\tRisc\n\tyy\tBisc)
#
#this func may also used to estimate gene expr, nexe version
#v1.0.0
#
#17:02
#2010/8/1

@filenames=@ARGV;

for $filename(@filenames)
{
	undef @buffer;#storing all the data
	undef @tag;#deciding how to put result, 0 normal, n first of n repeated  lines,-1 rest of the repeated lines


	undef $unique_number;
	undef $total_number;

	open F1, $filename;


	#read first column to memory
	while(<F1>)
	{
	@elem_1=split /\s+/,$_;
	unshift(@buffer, $elem_1[0]);
	}

	$total_number=@buffer;


	#search in column 1 buffer
	for($i=0; $i<$total_number; $i++)
	{
		if(!$tag[$i])
		{
			$tag[$i]=0;
			$unique_number++;
			for($j=$i+1;$j<$total_number;$j++)
			{
				if($buffer[$i] eq $buffer[$j])
				{
					$tag[$i]++;
					$tag[$j]=-1;
				}
			}
		}
	}


	print "unique counts: $unique_number\n";
	print "total counts: $total_number\n";
	
	$filename_out=$filename."_name_noredun";
	open FOUT, ">$filename_out"|| die "can't create output file\n";

	for($i=0; $i<$total_number; $i++)
	{
		for($j=0; $j <= $tag[$i] ; $j++)
		{
			print FOUT "$buffer[$i]\n";
		}
	}
	close FOUT;
}


