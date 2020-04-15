#!/usr/bin/perl -w
print <<"usage" if (@ARGV < 2);

                =======================
                ||                   ||
                ||  select_fasta.pl  ||
                ||                   ||
                ||      by: liuhui   ||
                =======================

FUNCTION:               select target section from fasta file

HOWTO:                  .pl limition_list fasta_file >result

Thu Mar 10 15:41:40 CST 2011

usage

exit if (@ARGV <2);


$list_name=$ARGV[0];
#$fasta_name=$ARGV[1];
my %list;
open LIST, $list_name or $list{$list_name}=1;
#open FASTA, $fasta_name or die;

if (fileno LIST)
{
while (<LIST>)
{
        chomp;
        s/transcript://;
        $list{(split)[0]} = 1 if($_ ne "");
}
}

$/='>';
open IN, $ARGV[1] or die;
while(<IN>)
{
    chomp;
    next if ($_ eq'');
    my ($name, $seq) = split/\n/, $_, 2;
#transcript:Zm00001d010101_T001
    $name=~s/\s.*//;
    $name=~s/_.*//;
    $hash{$name}= $seq if (exists $list{$name});
}

$/="\n";
for my $k(sort keys %hash)
{
    print ">$k\n$hash{$k}";
}
