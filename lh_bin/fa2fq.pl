#!/usr/bin/perl -w

$fname=$ARGV[0];

@e=split/\./, $fname;
$oname=$e[0].".fq";

open FIN, $fname;
open FOUT, ">$oname"||die "can't create output file\n";

while(<FIN>)
{
chomp;
$iden=substr $_, 0, 1;
if($iden && $iden eq ">")
{
$line=substr($_, 1);
chomp($seq=<FIN>);

print FOUT "\@$line\n$seq\n+\n";

for($i=0; $i<length($seq); $i++)
{
print FOUT "I";
}
print FOUT "\n";
}
}
