#!/usr/bin/perl -w

#only diff in the storage source with merge_redun.pl


#function: merged repeated lines by the first column
#
#usage: merge_redun.pl tobemergedfile(can be multi) >cout_file
#
#input format:
#1.tobemergedfile
#11 xx  Risc
#22 zz  fds
#11 yy  Bisc
#
#output specification:
#1.tobemergedfile_noredun
#11 xx  Risc
#   yy  Bisc
#22 zz  fds
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

  open F1, $filename;
  $filename_out=$filename."_noredun";
  open FOUT, ">$filename_out"|| die "can't create output file\n";

  undef $unique_number;
  undef $total_number;

  undef @repeated_namelist;

  undef @buffer;

  while(<F1>)
	{
	  push(@buffer, $_);
	}

	close F1;



#
#
#begin reading lines in buffer
#
#
  for $in1(@buffer)
  {
	undef $line_print;#final one

	$line_1=$in1;

	chomp($line_1);

	$total_number++;#current total number actually

	@elem_1=split /\s+/,$line_1;
	$line_1=join "\t", @elem_1;


	$name_1=shift(@elem_1);#1

	unless (&isfoundbefore())
	{
		$unique_number++;

		$line_print.="$line_1\n";
		$search_result=&search_repeated();
		if($search_result)
		{
			$line_print.="$search_result"; 
			unshift(@repeated_namelist, $name_1);
		}
		print FOUT $line_print;
	}

  }

  print "file: $filename\n";
  print "unique reads number: $unique_number\n";
  print"total reads number: $total_number\n";

}


#search global var name_1 in repeated_namelist, return result
sub isfoundbefore
{
	for $rep_name(@repeated_namelist)
	{
		return 1 if($rep_name eq $name_1);
	}
	0;
}


#global var used
#pre_limit: not part of a found repeated line
sub search_repeated
{
	my($result);

	for($i=$total_number; $i<@buffer; $i++)
	{
		@elem_2=split /\s+/,$buffer[$i];
		$name_2=shift(@elem_2);#1

		if($name_1 eq $name_2)
		{
			$result.= join "\t", ("",@elem_2,"\n");
		}
	}

	$result;
}

