#!/usr/bin/perl -w
#
#Author: Ruan Jue <ruanjue@genomics.org.cn>
#
use warnings;
use strict;

my $min_length = 0;
my $name = '';
my $seq = '';

while(<>){
        if(/^>(\S+)/){
                &print_scafftig($name, $seq) if($seq);
                $name = $1;
                $seq  = '';
        } else {
                chomp;
                $seq .= $_;
        }
}
&print_scafftig($name, $seq) if($seq);

1;

sub print_scafftig {
        my $name = shift;
        my $seq  = shift;
        my $id = 1;
        while($seq=~/([ATGCatgc]+)/g){
                my $s = $1;
                next if(length($s) < $min_length);
                print ">$name\_$id length=".length($s)."\n";
                while($s=~/(.{1,60})/g){
                        print "$1\n";
                }
                $id ++;
        }
}


