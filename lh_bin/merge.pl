#!/usr/bin/perl -w

#only work when current dir is where data located

#print @ARGV;automatically trans *.bz2 to filenames.bz2

#for($index=0; $index<@_; $index++)
#{
#  if(substr($ARGV[$index], 0, 1) eq "*")
  #{
 #   push(@files, glob($ARGV[$index]));
  #}
 # else 
 # {
 #   push(@files, $ARGV[$index]);
 # }
#}

@files=@ARGV;

#print @files;

foreach $filename(@files)
{
  open IN, $filename;
  print STDOUT "$filename\n-----------------------------------------------------------\n";
  while(<IN>)
  {
    print STDOUT;
  }
  print STDOUT "\n";
}
